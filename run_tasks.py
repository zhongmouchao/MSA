#!/usr/bin/env python3
"""
run_tasks.py -- Step 3 of the viroGym surveillance pipeline.

Turns task_plan.json into SLURM jobs, with the politeness-aware wave
scheduler:

  - one job per (reference group x FASTA file)
  - all files of one entry (group) are submitted together, in parallel
  - entries with MORE THAN ONE file run EXCLUSIVELY: no other entry is
    submitted until every job of that entry has left the queue
  - single-file entries are batched, up to --parallel-small entries at a
    time (default 10)
  - already-finished chunks (.msa_complete) and already-queued jobs are
    never resubmitted -- re-running this script is how you resume

Run from the pipeline root, AFTER plan_tasks.py:

    python3 run_tasks.py --dry-run          # show the wave plan, submit nothing
    python3 run_tasks.py                    # submit and babysit until done
    python3 run_tasks.py --no-wait          # submit wave 1 only is NOT what
                                            # this does: it submits every wave
                                            # back-to-back WITHOUT waiting;
                                            # only use if you accept many jobs
                                            # in the queue at once

Job mechanics (run_pipeline.py is NOT modified):
  - each job gets a staging dir  jobs/<group>/<chunk>/input/
    containing a symlink to the real FASTA plus the group config copied to
    <chunk>.txt  (this is how run_pipeline.py finds the right config)
  - outputs land in aligned/<group>/<chunk>/  with .msa_complete on success
  - job names: vg_<group>_<chunk>   (used for queue dedup and wave waiting)
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

# years that make no sense for recent-outbreak viruses; keeps garbage
# header dates out of the count tables
START_YEAR = {"SCV2": 2019}

JOB_PREFIX = "vg_"
MAX_JOB_NAME = 200  # safety; SLURM itself allows long names


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def queued_job_names(user):
    """Set of job names currently pending/running for this user, or None on error."""
    r = sh(["squeue", "-u", user, "-h", "-o", "%j"])
    if r.returncode != 0:
        return None
    return set(r.stdout.split())


def job_name(group_id, chunk_stem):
    name = f"{JOB_PREFIX}{group_id}_{chunk_stem}"
    return name[:MAX_JOB_NAME]


def chunk_stem(rel_path):
    return os.path.splitext(os.path.basename(rel_path))[0]


def stage_job(root, group, rel, fasta_abs):
    """Create jobs/<group>/<chunk>/input/ with fasta symlink + config copy."""
    gid = group["group_id"]
    stem = chunk_stem(rel)
    job_dir = os.path.join(root, "jobs", gid, stem)
    in_dir = os.path.join(job_dir, "input")
    os.makedirs(in_dir, exist_ok=True)

    link = os.path.join(in_dir, os.path.basename(rel))
    if os.path.islink(link) or os.path.exists(link):
        os.remove(link)
    os.symlink(fasta_abs, link)

    cfg_dst = os.path.join(in_dir, stem + ".txt")
    shutil.copyfile(group["config"], cfg_dst)
    return job_dir, in_dir


def write_manifest(root, group):
    """So scan_aligned.py can see this group from now on."""
    gdir = os.path.join(root, "aligned", group["group_id"])
    os.makedirs(gdir, exist_ok=True)
    manifest = {
        "group_id": group["group_id"],
        "ref_hash": group["group_id"].rsplit("_", 1)[-1],
        "ref_len": group["ref_len"],
        "config": group["config"],
        "motif_tags": [m["tag"] for m in group["motifs"]],
        "registered_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "source_folder_name": "run_tasks.py",
    }
    path = os.path.join(gdir, "alignment_manifest.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)


def chunk_done(root, gid, stem):
    return os.path.isfile(os.path.join(root, "aligned", gid, stem, ".msa_complete"))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", default="task_plan.json")
    ap.add_argument("--parallel-small", type=int, default=10,
                    help="single-file entries allowed in flight at once "
                         "(default 10; multi-file entries always run alone)")
    ap.add_argument("--poll", type=int, default=180,
                    help="seconds between queue checks while waiting for a "
                         "wave to finish (default 180)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the wave plan and submit nothing")
    ap.add_argument("--no-wait", action="store_true",
                    help="submit all waves without waiting between them "
                         "(impolite on shared clusters; for testing)")
    ap.add_argument("--only-group", default=None,
                    help="substring filter: run only matching group ids")
    args = ap.parse_args()

    root = os.getcwd()
    if not os.path.isfile(args.plan):
        print(f"[FATAL] {args.plan} not found -- run plan_tasks.py first")
        return 1
    with open(args.plan, encoding="utf-8") as fh:
        plan = json.load(fh)

    user = os.environ.get("USER") or os.environ.get("LOGNAME") or "unknown"

    # ---- collect pending jobs per group ------------------------------------
    groups_todo = []
    n_skip_done = 0
    for g in plan.get("groups", []):
        if args.only_group and args.only_group not in g["group_id"]:
            continue
        todo = []
        for f in g["fasta_files"]:
            stem = chunk_stem(f["rel_path"])
            if chunk_done(root, g["group_id"], stem):
                n_skip_done += 1
                continue
            todo.append(f)
        if todo:
            gg = dict(g)
            gg["files_todo"] = todo
            groups_todo.append(gg)

    if not groups_todo:
        print("nothing to do -- every planned chunk already has .msa_complete")
        return 0

    total_jobs = sum(len(g["files_todo"]) for g in groups_todo)
    print("=" * 72)
    print("RUN PLAN")
    print("=" * 72)
    print(f"groups with pending work : {len(groups_todo)}")
    print(f"jobs to submit           : {total_jobs}")
    print(f"chunks already finished  : {n_skip_done} (skipped)")
    print("-" * 72)

    # ---- wave construction --------------------------------------------------
    big = [g for g in groups_todo if len(g["files_todo"]) > 1]
    small = [g for g in groups_todo if len(g["files_todo"]) == 1]
    waves = [[g] for g in big]
    for i in range(0, len(small), args.parallel_small):
        waves.append(small[i:i + args.parallel_small])

    for wi, wave in enumerate(waves, 1):
        njobs = sum(len(g["files_todo"]) for g in wave)
        kinds = ",".join(f"{g['group_id']}({len(g['files_todo'])})" for g in wave)
        print(f"  wave {wi:>3}: {njobs:>4} job(s)  <- {kinds}")
    if args.dry_run:
        print("\n--dry-run: nothing submitted")
        return 0

    # ---- execute waves --------------------------------------------------------
    submitted_names = set()
    for wi, wave in enumerate(waves, 1):
        print()
        print("=" * 72)
        print(f"WAVE {wi}/{len(waves)}")
        print("=" * 72)

        # never resubmit something already in the queue
        queued = queued_job_names(user)
        if queued is None:
            print("[FATAL] squeue failed -- refusing to submit blindly. "
                  "Re-run when SLURM responds; finished chunks are skipped.")
            return 1

        wave_names = []
        for g in wave:
            write_manifest(root, g)
            binning = g.get("binning", "yearly")
            extra = []
            if binning == "monthly":
                extra.append("--monthly")
            sy = START_YEAR.get(g["virus_code"])
            if sy:
                extra += ["--start-year", str(sy)]

            for f in g["files_todo"]:
                gid = g["group_id"]
                stem = chunk_stem(f["rel_path"])
                name = job_name(gid, stem)
                if name in queued or name in submitted_names:
                    print(f"  SKIP (queued): {name}")
                    continue
                job_dir, in_dir = stage_job(root, g, f["rel_path"], f["path"])
                out_dir = os.path.join(root, "aligned", gid, stem)
                cmd = ["sbatch", "--job-name", name,
                       os.path.join(root, "vg_run_slurm.sh"),
                       "--input-dir", in_dir,
                       "--output-dir", out_dir,
                       "--dataset", os.path.basename(f["rel_path"])] + extra
                r = sh(cmd, cwd=root)
                if r.returncode == 0:
                    print(f"  submitted {name}  ({r.stdout.strip()})")
                    wave_names.append(name)
                    submitted_names.add(name)
                else:
                    print(f"  [WARNING] submission failed for {name}: "
                          f"{r.stderr.strip() or r.stdout.strip()}")
                time.sleep(1)  # gentle on the scheduler

        if args.no_wait:
            continue
        if not wave_names:
            print("  (nothing new submitted this wave; not waiting)")
            continue

        # wait until every job of this wave has left the queue
        print(f"  waiting for {len(wave_names)} job(s) to finish "
              f"(polling every {args.poll}s) ...")
        while True:
            time.sleep(args.poll)
            queued = queued_job_names(user)
            if queued is None:
                print("  squeue hiccup -- retrying next poll")
                continue
            remaining = [n for n in wave_names if n in queued]
            done_now = len(wave_names) - len(remaining)
            print(f"    {time.strftime('%H:%M:%S')}: {done_now}/"
                  f"{len(wave_names)} finished")
            if not remaining:
                break

    print()
    print("=" * 72)
    print("ALL WAVES SUBMITTED AND FINISHED (queue-wise)")
    print("=" * 72)
    # final honesty check against disk state
    missing = []
    for g in plan.get("groups", []):
        if args.only_group and args.only_group not in g["group_id"]:
            continue
        for f in g["fasta_files"]:
            stem = chunk_stem(f["rel_path"])
            if not chunk_done(root, g["group_id"], stem):
                missing.append(f"{g['group_id']} / {stem}")
    if missing:
        print(f"{len(missing)} chunk(s) have NO .msa_complete marker --")
        print("they failed or were cancelled. Inspect their logs, then simply")
        print("re-run this script; only missing chunks will be resubmitted:")
        for m in missing[:20]:
            print(f"  - {m}")
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more")
        return 1
    print("every planned chunk has .msa_complete -- alignment phase done")
    print("next: python3 aggregate_counts.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
