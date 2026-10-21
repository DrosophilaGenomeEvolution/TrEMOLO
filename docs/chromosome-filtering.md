# Chromosome-name filters

An [inventory of the matching source lines](chromosome-filter-usage.md) lists
configuration, validation and application sites, separating the current
workflow, historical workflow and independent utilities.

In the refactored `workflow/Snakefile`, `CHROM_KEEP` is a **Python regular
expression applied to the chromosome identifier alone**. It never sees the
FASTA description, a trailing space/tab, BED coordinates, TE labels or read
identifiers. The same matcher is used for Sniffles 1/2 VCF sequences/read names,
OUTSIDER BED merging and resident TE annotation.

For exact Drosophila main-chromosome names:

```yaml
PARAMS:
  OUTSIDER_VARIANT:
    TE_DETECTION:
      CHROM_KEEP: '^([23][RL]|[X4Y])$'
  TE_GENOME:
    CHROM_KEEP: '^([23][RL]|[X4Y])$'
```

These two settings have separate scopes: the first filters OUTSIDER calls,
the second filters resident TE annotation. Setting one does not silently set
the other. Their defaults remain `.` (keep all chromosome names).

The example accepts `2L`, `2R`, `3L`, `3R`, `X`, `4`, `Y` and excludes names
such as `2L_extra`, `X_random` and `scaffoldX`. A FASTA header `>2L description`
has identifier `2L`; the description is not searched. A name actually containing
an assembly suffix must be selected explicitly if it is wanted.

Matching uses `re.search`: add `^` and `$` for an exact name. Without `$`,
`^2L` intentionally accepts `2L_extra`. Legacy comma-separated lists such as
`^2L$,^X$` remain supported as alternatives. Commas inside character classes,
groups and quantifiers (for example `^chr[12]{1,2}$`) are not list separators.
All stages now use Python regex syntax rather than a mixture of Python and
`grep -E` syntax. BED rows are streamed and preserved unchanged after selecting
their first field.

Prefer YAML single quotes for regular expressions. Double quotes interpret
backslashes: `"\\s"` passes `\s`, whereas `"\s"` is invalid YAML. Regardless
of quoting, `\s` requires whitespace and will not match the end of a chromosome
identifier. Existing separator-based expressions should be replaced by name
expressions such as the example above; the pipeline does not rewrite them.

Invalid regexes fail before workflow jobs execute. Both Sniffles readers fail
explicitly if a VCF declares contigs but none matches the requested filter;
the error explains that only the name is searched and lists example names.
This avoids a misleading successful empty extraction with Sniffles 2. An empty
VCF whose contig names match still yields a valid empty result, and VCFs without
contig declarations are filtered on the chromosome field of their records.

The GFA export's separate `--chrom-regex` option already searches chromosome
names only. It accepts one Python expression (use `|` for alternatives), not
the legacy comma-list convention. See [GFA export](te-gfa.md).

The historical `run.snk`/`rules.snk` workflow is not the refactored entry point
and does not acquire the new BED-merging contract. Use `workflow/Snakefile`
for this behavior.
