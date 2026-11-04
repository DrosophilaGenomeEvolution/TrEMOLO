#!/usr/bin/env python3
"""Identify the source checkout, independently of the execution container."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def version_info(root=ROOT):
    root = Path(root).resolve()
    info = {"version": (root / "VERSION").read_text().strip(),
            "commit": None, "dirty": None}
    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(root)] + list(args), stderr=subprocess.DEVNULL,
            universal_newlines=True, timeout=5).strip()
    try:
        if Path(git("rev-parse", "--show-toplevel")).resolve() == root:
            info["commit"] = git("rev-parse", "HEAD")
            info["dirty"] = bool(git("status", "--porcelain", "--untracked-files=normal"))
    except (OSError, subprocess.SubprocessError):
        pass
    return info


def banner(info=None):
    info = version_info() if info is None else info
    commit = info["commit"] or "unavailable (source archive or Git unavailable)"
    state = " [modified checkout]" if info["dirty"] else ""
    return "TrEMOLO {} | commit {}{}".format(info["version"], commit, state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Machine-readable source identity")
    args = parser.parse_args()
    info = version_info()
    print(json.dumps(info, sort_keys=True) if args.json else banner(info))


if __name__ == "__main__":
    main()
