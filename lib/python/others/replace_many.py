#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def strip_prefix(value: str) -> str:
    for prefix in ("original=>", "original=", "pseudo=>", "pseudo="):
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def load_allowed_ids(path: str | None) -> set[str] | None:
    if path is None:
        return None

    allowed_ids = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value:
                allowed_ids.add(value)

    return allowed_ids


def load_mapping(term_file: str, allowed_ids: set[str] | None) -> dict[str, str]:
    mapping: dict[str, str] = {}

    with open(term_file, "r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            stripped = line.rstrip("\n")
            if not stripped:
                continue

            parts = stripped.split("\t")
            if len(parts) < 2:
                continue

            original = strip_prefix(parts[0])
            pseudo = strip_prefix(parts[1])

            if not pseudo or not original or pseudo == original:
                continue

            if allowed_ids is not None and pseudo not in allowed_ids:
                continue

            mapping[pseudo] = original.replace("\\", "_")

    return mapping


def replace_entire_file(target_file: Path, mapping: dict[str, str]) -> None:
    pattern = re.compile(
        r"\b(?:"
        + "|".join(sorted((re.escape(term) for term in mapping), key=len, reverse=True))
        + r")\b"
    )

    content = target_file.read_text(encoding="utf-8")
    updated = pattern.sub(lambda match: mapping[match.group(0)], content)

    if updated != content:
        target_file.write_text(updated, encoding="utf-8")


def replace_column(target_file: Path, mapping: dict[str, str], column: int) -> None:
    index = column - 1
    lines = target_file.read_text(encoding="utf-8").splitlines(keepends=True)
    updated_lines = []
    changed = False

    for line in lines:
        newline = "\n" if line.endswith("\n") else ""
        raw_line = line[:-1] if newline else line
        columns = raw_line.split("\t")

        if 0 <= index < len(columns) and columns[index] in mapping:
            columns[index] = mapping[columns[index]]
            changed = True

        updated_lines.append("\t".join(columns) + newline)

    if changed:
        target_file.write_text("".join(updated_lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replace multiple pseudo TE identifiers in one pass.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("term_file", help="Tab-separated file with original and pseudo TE names.")
    parser.add_argument("target_file", help="File to update in place.")
    parser.add_argument(
        "--ids-file",
        help="Optional file containing one pseudo TE identifier per line to restrict replacements.",
    )
    parser.add_argument(
        "--column",
        type=int,
        help="Replace only the exact value found in the given 1-based tab-separated column.",
    )

    args = parser.parse_args()

    allowed_ids = load_allowed_ids(args.ids_file)
    mapping = load_mapping(args.term_file, allowed_ids)
    if not mapping:
        return

    target_file = Path(args.target_file)
    if not target_file.exists() or target_file.stat().st_size == 0:
        return

    if args.column is None:
        replace_entire_file(target_file, mapping)
    else:
        replace_column(target_file, mapping, args.column)


if __name__ == "__main__":
    main()
