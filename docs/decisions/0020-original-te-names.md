# 0020 — Restore original TE names in public outputs

Status: accepted.

## Problem

Input preparation uses safe `TrEMOLOTE<N>` identifiers when a TE FASTA header
contains characters that conflict with downstream parsers. Variable TE outputs
and the migrated report did not restore the original labels. The legacy
`RENAME_PSEUDO_TE` rule performed repeated substitutions in existing results,
without declaring the modified files as outputs.

## Decision

`restore_te_names` reads `INPUT/te_headers.tsv` once and publishes declared
outputs using exact dictionary lookups in family columns. Each file is streamed
once. An identity mapping uses a byte-for-byte copy. There is no subprocess per
family and no global replacement in event identifiers, read identifiers or
chromosome names. Each destination is replaced atomically after a successful
export; source files remain unchanged.

The public `TE_INFOS.bed`, `TE_CALL_CANDIDATES.tsv` and root position BEDs retain
their existing paths. Their computational inputs are under `INTERNAL/`.
Detailed human-readable CSV/TSV/BED outputs are mirrored under `PUBLIC/`, keeping
their branch-relative paths: classified and combined calls, family counts,
frequencies, TSDs, insertion/deletion BEDs, integration audit and Liftoff audit.
Existing detailed branch files remain computational outputs with prepared
identifiers. `PUBLIC/name-restoration.json` lists each source, destination and
number of restored cells. Only enabled branches are exported. Resident TE
annotation already restores names and is not passed through a second mapping.

Original labels containing spaces, quotes, `|` or `:` are preserved using TSV
quoting. Read the tables with a TSV CSV parser. In public `TE|ID` cells, the final
`|` separates the event identifier from the original family name. The report
and insertion structure module use this convention. The structure module also
maps raw BLAST subject identifiers when the public call table's sibling
`INPUT/te_headers.tsv` is available, preserving links to reported families.

Input preparation now treats whitespace and parser-sensitive punctuation as
requiring pseudonyms. Tabs in names are rejected explicitly. The mapping itself
uses TSV quoting so quotes in original FASTA headers round-trip correctly.
Report data schema is 1.5.0. The restoration benchmark is included in the report's
workflow timing chart.

## Validation

Unit coverage checks exact lookups (`TrEMOLOTE1` versus `TrEMOLOTE10`), original
names colliding with prepared identifiers, compound names, quote round trips,
unchanged provenance identifiers, repeatability, source immutability, identity
byte preservation, empty inputs and atomic failure behavior. Report and module
checks verify that a family containing `|` retains its original name and event
identifier. Full container workflows exercise ordinary and unsafe TE banks,
and Sniffles2 with clipped-read detection disabled.

Container validation passed 143 workflow tests, 8 module tests and 6 locus
model tests. Two fresh 71-job runs produced 529 calls each; a fresh Sniffles2
run with clipped reads disabled produced 19 calls in 43 jobs. Integration and
report regression checks passed. The ordinary bank's public call tables and
reconstructed genomes remain byte-identical to the previous migrated outputs.
Changing only FASTA names preserved coordinates, frequencies and both genomes;
public scientific columns match their internal sources exactly. On this fixture,
restoring 3,119 cells across 33 files took about 0.12 seconds.

The safe-versus-unsafe bank comparison also exposed an existing TSD preparation
limitation: `prepare_outsider_tsd.py` deliberately preserves the historical
`matched_direct[1:]` behavior, dropping the first matched direct candidate.
Changing family names changes classification ordering, so `TrEMOLO.INS.84`
received `NONE` for its TSD and TE-size fallback with the pseudonymized bank.
This happens before name restoration. Removing that compatibility behavior is
a separate scientific correction; this change preserves the internal result.
