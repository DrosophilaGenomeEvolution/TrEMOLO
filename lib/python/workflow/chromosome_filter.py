"""One Python-regex contract for chromosome identifiers in VCF, BED and BLAST."""
import argparse
from functools import lru_cache
from pathlib import Path
import re
import sys


def split_patterns(expression):
    """Keep legacy comma-separated patterns without splitting regex quantifiers."""
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("CHROM_KEEP must be a nonempty regular expression")
    patterns, start = [], 0
    escaped = bracket = False
    parentheses = braces = 0
    for index, character in enumerate(expression):
        if escaped:
            escaped = False
            continue
        if character == "\\":
            escaped = True
        elif bracket:
            if character == "]":
                bracket = False
        elif character == "[":
            bracket = True
        elif character == "(":
            parentheses += 1
        elif character == ")":
            parentheses = max(0, parentheses - 1)
        elif character == "{":
            braces += 1
        elif character == "}":
            braces = max(0, braces - 1)
        elif character == "," and not parentheses and not braces:
            patterns.append(expression[start:index].strip())
            start = index + 1
    patterns.append(expression[start:].strip())
    if any(not pattern for pattern in patterns):
        raise ValueError("CHROM_KEEP contains an empty comma-separated pattern")
    return patterns


class ChromosomeFilter:
    def __init__(self, expression):
        self.patterns = tuple(re.compile(pattern) for pattern in split_patterns(expression))

    def search(self, chromosome):
        for pattern in self.patterns:
            match = pattern.search(chromosome)
            if match is not None:
                return match
        return None


def compile_chromosome_filter(expression):
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("CHROM_KEEP must be a nonempty regular expression")
    return _cached_filter(expression)


@lru_cache(maxsize=64)
def _cached_filter(expression):
    try:
        return ChromosomeFilter(expression)
    except re.error as error:
        raise ValueError(f"Invalid CHROM_KEEP regular expression {expression!r}: {error}") from error


def validate_vcf_chromosomes(path, pattern):
    names = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            match = re.match(r"##contig=<ID=([^,>]+)", line)
            if match:
                names.append(match.group(1))
    if names and not any(pattern.search(name) for name in names):
        examples = ", ".join(names[:12])
        raise ValueError(
            "CHROM_KEEP matches no chromosome name in the VCF. "
            "The regex is applied to the name alone, without whitespace, description or other columns. "
            f"Available names: {examples}" + (f" … ({len(names)} contigs)" if len(names) > 12 else "")
        )


def filter_rows(rows, output, pattern):
    for line in rows:
        fields = line.split(None, 1)
        if fields and pattern.search(fields[0]):
            output.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regex", required=True, help="Python regex applied to column 1 only")
    parser.add_argument("--input", type=Path, help="BED/TSV file; stdin by default")
    args = parser.parse_args()
    try:
        pattern = compile_chromosome_filter(args.regex)
    except ValueError as error:
        parser.error(str(error))
    if args.input:
        with args.input.open() as handle:
            filter_rows(handle, sys.stdout, pattern)
    else:
        filter_rows(sys.stdin, sys.stdout, pattern)


if __name__ == "__main__":
    main()
