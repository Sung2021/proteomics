#!/usr/bin/env python3
"""
run_pipeline.py — Proteomics pipeline orchestrator
Mirrors the pattern used in lung-tme-deconv-profiler.

Usage:
  python python/run_pipeline.py                          # uses mode from config.yaml
  python python/run_pipeline.py --mode group             # Branch A only
  python python/run_pipeline.py --mode dose_response     # Branch B (05b + 06)
  python python/run_pipeline.py --steps 1 2 3            # specific steps
  python python/run_pipeline.py --rscript /path/to/Rscript
"""

import argparse
import subprocess
import sys
import os
import yaml


# --------------------------------------------------------------------------- #
# Step registry
# --------------------------------------------------------------------------- #
COMMON_STEPS = [
    (1, "R/01_load_spectra.R"),
    (2, "R/02_preprocess.R"),
    (3, "R/03_qc.R"),
    (4, "R/04_feature_extraction.R"),
]

BRANCH_A = [
    (5, "R/05a_limpa_dea.R"),
]

BRANCH_B = [
    (5, "R/05b_limpa_prefilter.R"),
    (6, "R/06_dromics_bmd.R"),
]


def load_mode_from_config(config_path: str = "config.yaml") -> str:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    return cfg.get("study", {}).get("mode", "group")


def build_steps(mode: str, step_filter: list[int] | None) -> list[tuple[int, str]]:
    branch = BRANCH_A if mode == "group" else BRANCH_B
    all_steps = COMMON_STEPS + branch
    if step_filter:
        all_steps = [(n, s) for n, s in all_steps if n in step_filter]
    return all_steps


def run_step(step_num: int, script: str, rscript: str) -> bool:
    print(f"\n{'='*60}")
    print(f"  Step {step_num}: {script}")
    print(f"{'='*60}")

    result = subprocess.run([rscript, "--vanilla", script], capture_output=False)

    if result.returncode != 0:
        print(f"\n[ERROR] Step {step_num} failed (exit code {result.returncode})")
        print(f"Re-run manually: {rscript} --vanilla {script}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="Proteomics pipeline orchestrator")
    parser.add_argument(
        "--mode", choices=["group", "dose_response"],
        help="Override study.mode in config.yaml"
    )
    parser.add_argument(
        "--steps", nargs="+", type=int, metavar="N",
        help="Run only these step numbers (e.g. --steps 1 2 3)"
    )
    parser.add_argument(
        "--rscript", default="Rscript",
        help="Path to Rscript executable (default: Rscript)"
    )
    args = parser.parse_args()

    # Resolve working directory to repo root
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)

    mode = args.mode or load_mode_from_config()
    steps = build_steps(mode, args.steps)

    print(f"\nProteomics pipeline")
    print(f"Mode   : {mode}")
    print(f"Steps  : {[n for n, _ in steps]}")

    for step_num, script in steps:
        ok = run_step(step_num, script, args.rscript)
        if not ok:
            sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  Pipeline complete ({mode} mode)")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
