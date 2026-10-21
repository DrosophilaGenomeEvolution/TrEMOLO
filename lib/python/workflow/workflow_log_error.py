"""Select one bounded diagnostic from fresh tool logs after a job fails."""
import argparse
from pathlib import Path
import re


def diagnostic(paths, since=0):
    for path in reversed(paths):
        path = Path(path)
        try:
            if path.stat().st_mtime < since:
                continue
            with path.open("rb") as handle:
                handle.seek(0, 2)
                size = handle.tell()
                handle.seek(max(0, size - 32768))
                lines = handle.read(32768).decode("utf-8", errors="replace").splitlines()
            if size > 32768:
                lines = lines[1:]  # The first line may be an incomplete chunk.
        except OSError:
            continue
        lines = [re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", line).strip() for line in lines]
        candidates = [line for line in lines if re.search(
            r"\b(?:[A-Za-z_]*(?:Error|Exception)|ERROR|FATAL)\b|\berror:|\bfailed\b|not found|No such file",
            line, re.IGNORECASE) and not line.startswith("Traceback")]
        if candidates:
            return " ".join(candidates[-1].split())[:400]
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--since", type=float, default=0)
    parser.add_argument("logs", nargs="*")
    args = parser.parse_args()
    message = diagnostic(args.logs, args.since)
    if message:
        print(message)


if __name__ == "__main__":
    main()
