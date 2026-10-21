# Workflow progress messages

The refactored `workflow/Snakefile` emits short messages for major scientific
jobs and report generation. Index-only jobs retain Snakemake's usual display.
Messages go to stderr (the job's inherited console/scheduler stream) and are
appended to `WORK_DIRECTORY/log/progress.log`. They do not enter FASTA, BED,
JSON or other scientific stdout streams.

Each line includes a timestamp with timezone offset, the work-directory name,
the rule name and a status:

```text
… [TrEMOLO][analysis][extract_outsider_sniffles_sequences][START] OUTSIDER: extracting Sniffles insertion sequences
… [TrEMOLO][analysis][extract_outsider_sniffles_sequences][STEP] Reading sniffles1 VCF and applying chromosome-name and sequence filters
… [TrEMOLO][analysis][extract_outsider_sniffles_sequences][STEP] Recovering read insertions at VCF events and recording their provenance
… [TrEMOLO][analysis][extract_outsider_sniffles_sequences][DONE] OUTSIDER: extracting Sniffles insertion sequences
```

`START` and `DONE` bracket a major job. A few compound jobs also emit `STEP`
at phase boundaries: BAM preparation, Sniffles sequence extraction, hard-clip
preparation, INSIDER insertion/deletion BLAST and OUTSIDER evidence merging.
TE_GENOME target names are included in their descriptions. There is no message
per read, HSP, chromosome block or worker.

On a nonzero exit, `ERROR` identifies the last phase and lists the detailed
log path. When a diagnostic can be recognized in a fresh log, a second `ERROR`
line contains at most 400 characters from that diagnostic. Only the last 32 KiB
of each tool log are inspected, and only on failure. The original command's
exit status is preserved; an error does not produce a misleading `DONE`.
Unrecognized errors still have the phase, exit status and log path; consult
the detailed log for the full cause. No heuristic can reliably summarize every
external tool's error format.

Messages from parallel jobs may interleave. Their labels identify the source;
they are not a single sequential progress bar or a percentage of elapsed time.
In cluster execution, console messages follow the scheduler's normal stderr
routing. `progress.log` requires the usual shared work-directory access and
contains only jobs using the updated logging helper.

The file is append-only across resumed runs. Existing successful jobs are not
retroactively logged. Fine-grained output remains in each existing per-rule
log. Snakemake's own scheduling/error messages remain available.

Normal messages use Bash builtins, with no read/BAM/FASTA rescan and no timer
or polling loop. The diagnostic Python helper runs only after failure. No
runtime overhead percentage has been measured on a complete real dataset.

Implementation: `lib/bash/workflow_logging.sh` and
`lib/python/workflow/workflow_log_error.py`.
