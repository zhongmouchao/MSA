#!/usr/bin/env python3
"""
main_run.py -- one-command entry point for the viroGym surveillance pipeline.

Runs the steps in order, stopping at the first failure:

    check_registry.py    validate Virogym_Benchmark/ against the registry CSV
    scan_inputs.py       inventory RAW_FASTA_DATA/
    scan_aligned.py      inventory aligned/ (resume state)
    plan_tasks.py        build task_plan.json + planned_configs/
    run_tasks.py         submit SLURM waves (skipped with --no-submit)
    aggregate_counts.py  merge counts -> training_data/ (skipped with
                         --no-aggregate; harmless when counts don't exist
                         yet -- those motifs are reported as 'no count files')

Typical use from the pipeline root:

    python3 main_run.py --dry-run           # validate + plan, show waves
    python3 main_run.py                     # full run: submit + wait + aggregate
    python3 main_run.py --no-submit         # validate + plan + aggregate only
                                            # (e.g. day 2 while jobs still run)

Extra flags after '--' are passed to run_tasks.py, e.g.:

    python3 main_run.py -- --poll 300 --parallel-small 8
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run_step(script, extra=None, cwd=None):
    print()
    print("#" * 72)
    print(f"# {script} {' '.join(extra or [])}")
    print("#" * 72)
    r = subprocess.run([sys.executable, os.path.join(HERE, script)] +
                       (extra or []), cwd=cwd or os.getcwd())
    if r.returncode != 0:
        print(f"\n[STOP] {script} exited with code {r.returncode}. "
              f"Fix the problem and re-run main_run.py -- "
              f"finished steps are skipped automatically.")
        return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="plan and show waves, submit nothing")
    ap.add_argument("--no-submit", action="store_true",
                    help="skip the SLURM submission step entirely")
    ap.add_argument("--no-aggregate", action="store_true",
                    help="skip the aggregation step")
    args, passthrough = ap.parse_known_args()
    if passthrough and passthrough[0] == "--":
        passthrough = passthrough[1:]

    steps = [
        ("check_registry.py", []),
        ("scan_inputs.py", []),
        ("scan_aligned.py", []),
        ("plan_tasks.py", []),
    ]
    for script, extra in steps:
        if not run_step(script, extra):
            return 1

    if not args.no_submit:
        extra = list(passthrough)
        if args.dry_run:
            extra.append("--dry-run")
        if not run_step("run_tasks.py", extra):
            return 1

    if not args.no_aggregate:
        if not run_step("aggregate_counts.py", []):
            return 1

    print()
    print("#" * 72)
    print("# PIPELINE PASS COMPLETE")
    print("#" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
