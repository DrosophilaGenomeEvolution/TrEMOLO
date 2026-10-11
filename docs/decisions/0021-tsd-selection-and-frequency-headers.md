# 0021 — Keep all TSD candidates and document frequency tables

Status: accepted. Supersedes the candidate-selection compatibility behavior in
0003 and the headerless frequency-table layout in 0004.

## TSD correction

TSD preparation previously dropped the first retained direct-insertion query,
then the first combined Sniffles/direct candidate, and discarded the resulting
list if it contained at most one row. These were historical shell header-removal
artifacts, not biological filters. They could make TSD recovery depend on TE
family names and could suppress TSD detection for a single insertion.

Every selected candidate is now retained. Direct candidates require exact event
membership in the merged-call BED; Sniffles candidates require exact query
membership in their filtered call table. Membership sets are separate per source,
so a query selected for one source cannot select an unretained query in another.
Dictionary/set lookups replace repeated regular-expression scans of entire rows.
Flank eligibility checks, the existing flank convention and motif algorithm
continue to determine which candidates can produce a TSD.

## Frequency contract

Both OUTSIDER frequency TSVs have an identical eleven-column header. Their
column order and numeric body are unchanged. Empty results contain the header
only. `build_te_infos.index_frequency` recognizes the exact schema header and
skips it; legacy headerless tables remain accepted. Public name restoration
uses the named `te_family` column. Frequency regression comparisons remove only
the recognized header before comparing data bytes and counting events.

The column descriptions, INS/DEL semantics, percentage units, support corrections,
and clipped-read behavior are documented in `docs/frequency-tables.md`.
The workflow copies that document to `OUTSIDER/FREQUENCY/README.md` and
`PUBLIC/OUTSIDER/FREQUENCY/README.md`, next to the corresponding tables.
`SV_SIZE.tsv` is a computational size table and retains its existing format.

## Validation

Tests cover the first candidate, a singleton in each source, exact event matching,
source isolation, deterministic flank preparation, header-only outputs and
byte-identical `TE_INFOS` generation from headerless versus headed frequency
tables. Complete safe/unsafe-bank and Sniffles2 runs validate scientific data,
public names, integration and report consumers together.

Validation passed 146 workflow tests, 8 module tests and 6 locus tests. Two
71-job runs and a 43-job Sniffles2 run completed from fresh outputs. Both
frequency tables' data bytes match the preceding run exactly after removing
the new header, in all three configurations. Both reconstructed genomes also
remain byte-identical. The documentation is present beside source and public
tables. Report and integration regression checks pass.

TSD preparation retains 28 candidates instead of 26 on the full fixture.
With unsafe family names, `TrEMOLO.INS.84` recovers `GTCC` and its TE length,
raising TSD calls from 8 to 9. All 529 final calls agree between ordinary and
unsafe banks, including TSD and TE-size fields, after normalizing labels and
family-dependent legacy TrEMOLO identifiers.
