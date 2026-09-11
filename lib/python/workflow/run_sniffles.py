#!/usr/bin/env python3
"""Run a selected Sniffles generation without changing the legacy command."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def command(caller, executable, version, bam, reference, vcf, threads, support):
    if caller == "sniffles1":
        if version not in ("1.0.10", "1.0.11", "1.0.12", "1.0.12b"):
            raise ValueError("Expected Sniffles 1.0.10–1.0.12b, found " + version)
        report = "--report_seq" if version == "1.0.10" else "--report-seq"
        return [executable, "-t", str(threads), report, "-s", str(support),
                "-m", str(bam), "-v", str(vcf), "-n", "-1"]
    if not version.startswith("2."):
        raise ValueError("Expected Sniffles 2, found " + version)
    return [executable, "--input", str(bam), "--vcf", str(vcf),
            "--reference", str(reference), "--threads", str(threads),
            "--minsupport", str(support), "--output-rnames", "--allow-overwrite"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", required=True, choices=("sniffles1", "sniffles2"))
    parser.add_argument("--executable", required=True)
    for name in ("bam", "reference", "vcf", "metadata"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--min-support", type=int, default=1)
    args = parser.parse_args()
    if args.min_support < 1 or args.threads < 1:
        parser.error("threads and min-support must be positive integers")
    probe = subprocess.run([args.executable, "-h" if args.caller == "sniffles1" else "--version"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    match = re.search(r"(?:Version\s*:?\s*|Sniffles2?\s*v?)(\d+\.\d+(?:\.\d+)?[a-z]?)",
                      probe.stdout, re.IGNORECASE)
    if not match:
        raise ValueError("Cannot identify Sniffles version: " + probe.stdout)
    version = match.group(1)
    argv = command(args.caller, args.executable, version, args.bam, args.reference,
                   args.vcf, args.threads, args.min_support)
    args.vcf.parent.mkdir(parents=True, exist_ok=True)
    print("Running " + repr(argv), flush=True)
    subprocess.run(argv, check=True)
    if not args.vcf.is_file() or not args.vcf.stat().st_size:
        raise ValueError("Sniffles did not produce a VCF header")
    args.metadata.write_text(json.dumps({"caller": args.caller, "version": version,
                                         "command": argv}, indent=2) + "\n")


if __name__ == "__main__":
    main()
