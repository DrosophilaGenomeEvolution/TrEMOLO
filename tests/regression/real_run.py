#!/usr/bin/env python3
"""Capture/check real-run outputs and record tools inside the execution environment."""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN_TABLES = ("TE_INFOS.bed", "TE_CALL_CANDIDATES.tsv", "TE_CALL_CANDIDATE_EVIDENCE.tsv")
GRAPH_FILES = ("te-graph.gfa", "sequence-alternatives.tsv", "loci.tsv", "excluded-events.tsv", "walks.tsv")
GRAPH_SETTINGS = ("reference_name", "construction", "locus_window", "thresholds", "selected_tes", "chrom_regex", "only_called_chromosomes")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def signature(path, tabular=False):
    result = dict(sha256=sha256(path))
    if tabular:
        hashes = []
        with Path(path).open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None:
                raise ValueError(f"Missing table header: {path}")
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f"Invalid row width: {path}")
                hashes.append(hashlib.sha256(json.dumps(row, ensure_ascii=False, separators=(",", ":")).encode()).digest())
        result.update(header=header, rows=len(hashes), semantic_sha256=hashlib.sha256(b"".join(sorted(hashes))).hexdigest())
    return result


def named_path(value):
    name, separator, path = value.partition("=")
    if not separator or not name or not path:
        raise argparse.ArgumentTypeError("Use NAME=/path")
    return name, Path(path)


def unique_paths(values):
    paths = {}
    for name, path in values:
        if name in paths:
            raise ValueError("Duplicate name: " + name)
        paths[name] = path
    return paths


def capture(run=None, graphs=(), config=None, versions=None, inputs=(), command=None):
    result = dict(schema_version="1.0.0", run=None, graphs={}, inputs={}, configuration=None,
                  tool_versions=None, command=command)
    if run is not None:
        run = Path(run)
        if not (run / "TE_INFOS.bed").is_file():
            raise FileNotFoundError(run / "TE_INFOS.bed")
        result["run"] = dict(source=str(run.resolve()), outputs={name: signature(run / name, True)
            for name in RUN_TABLES if (run / name).is_file()})
    for name, directory in unique_paths(graphs).items():
        manifest = json.loads((directory / "manifest.json").read_text())
        result["graphs"][name] = dict(source=str(directory.resolve()),
            settings={key: manifest.get(key, False if key == "only_called_chromosomes" else None) for key in GRAPH_SETTINGS},
            recorded_inputs=manifest.get("inputs", []),
            outputs={filename: signature(directory / filename, filename.endswith(".tsv")) for filename in GRAPH_FILES})
    for name, path in unique_paths(inputs).items():
        result["inputs"][name] = dict(source=str(path.resolve()), sha256=sha256(path))
    if config is not None:
        result["configuration"] = dict(source=str(config.resolve()), sha256=sha256(config), text=config.read_text())
    if versions is not None:
        recorded = json.loads(versions.read_text())
        tools = recorded.get("tools", {})
        if not isinstance(tools, dict) or not tools or any(
                not isinstance(value, dict) or value.get("returncode") != 0 for value in tools.values()):
            raise ValueError("Version record must contain successfully queried tools")
        result["tool_versions"] = dict(source=str(versions.resolve()), data=recorded)
    if result["run"] is None and not result["graphs"]:
        raise ValueError("Supply --run and/or --gfa NAME=DIR")
    missing = []
    for field in ("run", "configuration", "tool_versions", "inputs", "command"):
        if not result[field]:
            missing.append(field)
    result["reproducibility"] = dict(complete=not missing, missing=missing,
        scope="full_run_record" if not missing else "available_outputs_only")
    return result


def differences(expected, actual, prefix="", strict=False):
    """Paths and capture provenance are descriptive, never comparison keys."""
    if prefix == "tool_versions.data.pipeline":
        # A regression run intentionally exercises a new pipeline revision.
        return []
    changes = []
    ignored = {"source", "recorded_inputs", "reproducibility", "command"}
    if isinstance(expected, dict) and isinstance(actual, dict):
        for key in sorted(set(expected) | set(actual)):
            if key in ignored or (key == "sha256" and "semantic_sha256" in expected and not strict):
                continue
            label = f"{prefix}.{key}" if prefix else key
            if key not in actual or key not in expected:
                changes.append(label + ": missing or unexpected field")
            else:
                changes.extend(differences(expected[key], actual[key], label, strict))
    elif expected != actual:
        changes.append(prefix + ": changed")
    return changes


def record_versions(tools):
    result = dict(python=platform.python_version(), tools={})
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, timeout=10)
        status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, timeout=10)
        result["pipeline"] = dict(git_commit=revision.stdout.strip() if revision.returncode == 0 else None,
                                  dirty=bool(status.stdout.strip()) if status.returncode == 0 else None)
    except (OSError, subprocess.TimeoutExpired) as error:
        result["pipeline"] = dict(git_commit=None, dirty=None, error=str(error))
    for name, executable in unique_paths(tools).items():
        flag = "-version" if name == "blastn" else "--version"
        command = [str(executable), flag]
        try:
            process = subprocess.run(command, capture_output=True, text=True, timeout=20)
            result["tools"][name] = dict(command=command, returncode=process.returncode,
                                        stdout=process.stdout.strip(), stderr=process.stderr.strip())
        except (OSError, subprocess.TimeoutExpired) as error:
            result["tools"][name] = dict(command=command, returncode=None, error=str(error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("snapshot", "check"):
        child = sub.add_parser(action)
        child.add_argument("--run", type=Path)
        child.add_argument("--gfa", action="append", default=[], type=named_path)
        child.add_argument("--config", type=Path)
        child.add_argument("--versions", type=Path)
        child.add_argument("--input", action="append", default=[], type=named_path)
        child.add_argument("--command", help="Exact pipeline command, recorded without executing it")
        if action == "snapshot":
            child.add_argument("--output", type=Path, required=True)
        else:
            child.add_argument("--baseline", type=Path, required=True)
            child.add_argument("--strict", action="store_true", help="Also require byte-identical tables")
    version_parser = sub.add_parser("versions")
    version_parser.add_argument("--tool", action="append", required=True, type=named_path,
                                help="NAME=executable; run inside the same container/job as the analysis")
    version_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "versions":
        data = record_versions(args.tool)
        args.output.write_text(json.dumps(data, indent=2) + "\n")
        return int(any(t["returncode"] != 0 for t in data["tools"].values()))
    if args.action == "check":
        expected = json.loads(args.baseline.read_text())
        if args.run is None and expected["run"]:
            args.run = Path(expected["run"]["source"])
        if not args.gfa:
            args.gfa = [(n, Path(g["source"])) for n, g in expected["graphs"].items()]
        if not args.input:
            args.input = [(n, Path(i["source"])) for n, i in expected["inputs"].items()]
        if args.config is None and expected["configuration"]:
            args.config = Path(expected["configuration"]["source"])
        if args.versions is None and expected["tool_versions"]:
            args.versions = Path(expected["tool_versions"]["source"])
    data = capture(args.run, args.gfa, args.config, args.versions, args.input, args.command)
    if args.action == "snapshot":
        if args.output.exists():
            parser.error("Baseline already exists; use a new filename after reviewing the change")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
        print("Baseline recorded. Reproducibility:", json.dumps(data["reproducibility"]))
        return 0
    changes = differences(expected, data, strict=args.strict)
    for message in changes:
        print(message)
    print("Regression check FAILED" if changes else "Regression check passed")
    return int(bool(changes))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
