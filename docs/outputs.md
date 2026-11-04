# Output files

Outputs depend on enabled branches. Start with these files:

| Path under WORK_DIRECTORY | Purpose |
| --- | --- |
| `TE_INFOS.bed` | Final variant-associated calls; 15-column legacy-compatible schema |
| `TE_CALL_CANDIDATES.tsv` | Alternative TE annotations for ambiguous reported calls |
| `TE_CALL_CANDIDATE_EVIDENCE.tsv` | Candidate sequence/read evidence and provenance |
| `REPORT/report.html` | Self-contained main report |
| `PUBLIC/` | Detailed exports with original TE names; internal calculation files retain prepared IDs |
| `input_manifest.json`, `INPUT/te_headers.tsv` | Input checksums and original TE header mapping |
| `log/progress.log`, `log/run-version.json` | Workflow progress and source identity for the latest execution |
| `TE_GENOME/GENOME/`, `TE_GENOME/REFERENCE/` | Optional resident-copy tables, BED and GFF3 |
| `UNMAPPED_TE/` | Optional TE evidence in reads without genomic alignment |
| `OUTSIDER/FREQUENCY/` | Read support and frequency tables |
| `OUTSIDER/TE_TOWARD_GENOME/` | Optional reconstructed genomes and integration audit |
| `OUTSIDER/INSIDER_VR/LIFT_OFF_AUDIT.tsv` | Optional per-insertion reference projection audit |

Module reports and GFA graphs are explicit post-processing outputs; `all` does
not build them. See [modules](../README.md#modules) and [GFA](te-gfa.md).
For TSV frequency columns, use [frequency tables](frequency-tables.md).
For biological interpretation, see [interpreting results](interpreting-results.md).

## Main outputs

The most useful output files are :

* The html report in **your_work_directory/REPORT/report.html** with interactive summaries and call tables

The output file **your_work_directory/TE_INFOS.bed** gathers all the necessary information.

|      chrom      |  start   | end      |   TE\|ID   |   strand  |    TSD     | pident | psize_TE | SIZE_TE |      NEW_POS     |  FREQ (%) | FREQ_WITH_CLIPPED (%) | SV_SIZE | ID_TrEMOLO  | TYPE |
| --------------- | -------- | -------- | ---------- | --------- | ---------- | --------- | -------- | ------- | ---------------- | ------- | -------------- | -------------- | -------------- | -------------- |
|  2R_RaGOO_RaGOO | 16943971 | 16943972 | roo\|svim.INS.175  |     +     |   GTACA   | 97.026 | 99.2  | 9006    | 16943978 | 28.5714 |    28.5714     | 9000  |  TE_ID_OUTSIDER.94047.INS.107508.0  | INS |
|  X_RaGOO_RaGOO  | 21629415 | 21629416 | ZAM\|Assemblytics_w_534  |     -     |   CGCG  | 98.6 | 90.5  | 8435    | 21629413         | 11.1111 |    10.0000   | 8000  | TE_ID_INSIDER.77237.Repeat_expansion.8 | Repeat_expansion |


 1.    `chrom` : chromosome
 2.    `start` : start position for the TE
 3.    `end` : end position for the TE
 4.    `TE|ID` :   TE name and ID in **SV.vcf**,**SV_SOFT.vcf**,**HARD.fasta** and **SV_INS_CLUST.bed** (for OUTSIDER) or **assemblytics_out.Assemblytics_structural_variants.bed** (for INSIDER)
 5.    `strand` :  strand of the TE
 6.    `TSD` : TSD SEQUENCE
 7.    `pident` : percentage of identical matches with TE
 8.    `psize_TE` : percentage of size with TE in database
 9.    `SIZE_TE` :  TE size
 10.   `NEW_POS` :  position corrected with calculated TSD (only for OUTSIDER)
 11.   `FREQ`  : frequency, normalized
 12.   `FREQ_WITH_CLIPPED`  : frequency with clipped read (OUTSIDER only)
 13.   `SV_SIZE`  : size of the structural variant (may be larger than the size of the TE)
 14.   `ID_TrEMOLO`  : TrEMOLO ID of the TE
 15.   `TYPE`  : type of insertion can be HARD,SOFT (Warning : HARD, SOFT are often false positives),INS,INS_DEL... (INS_DEL is an insertion located on a deletion of the assembly)

