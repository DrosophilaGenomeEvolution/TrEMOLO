# ADR 0001: Pluggable INSIDER structural-variant backends

- Status: accepted for incremental implementation
- Date: 2026-07-16

## Context

TrEMOLO currently embeds a customized Assemblytics workflow. It is reproducible
and its output is required for backward compatibility, but it should not remain
the only possible interpretation of assembly-to-reference differences.

SyRI provides a broader synteny and rearrangement model for chromosome-level
assemblies. Other assembly callers may also provide useful independent evidence.
Changing the default caller immediately would invalidate established TrEMOLO
results before biological concordance has been measured.

## Decision

1. Keep `assemblytics` as the default and legacy-compatible INSIDER backend.
2. Design the workflow around a selectable `SV_CALLER` backend.
3. Reserve `syri` as the next experimental backend.
4. Normalize every backend into a common internal SV table before TE detection.
5. Allow a future `consensus` mode that records caller agreement instead of
   silently selecting one result.
6. Attach every normalized SV and locus assignment to its caller evidence.

Planned configuration:

```yaml
PARAMS:
  INSIDER_VARIANT:
    SV_CALLER: assemblytics  # assemblytics, syri, consensus
```

Only `assemblytics` is implemented now. Selecting an unavailable backend must
fail during configuration validation, never fall back silently.

## Consequences

- The current regression oracle remains valid.
- TE detection rules must consume normalized SVs rather than caller-specific
  files once a second backend is introduced.
- The locus model can retain multiple pieces of caller evidence for one event.
- SyRI adoption requires dedicated chromosome-level and rearrangement test data.

## Deferred validation

Before changing the default, compare recall, breakpoint agreement, TE-family
assignment, complex-event representation, and false merges on synthetic and
manually reviewed loci.
