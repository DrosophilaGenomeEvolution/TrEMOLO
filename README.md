

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5391)


<img src="images/TrEMOLO9.png">

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
    - [Using Git](#git)
    - [Using Singularity](#singularity)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output files](#output)
- [Citation & Licence](#citation)


# TrEMOLO<a name="introduction"></a>
**Transposable Elements MOvement detection using LOng reads**

TrEMOLO used long reads, raw or their assemblies to detect
- Global TE variations between two assembled genomes
- Populational/somatic variation in TE insertions



## Requirements<a name="requirements"></a>

- For both approaches
  - Python 3.6+
- For Global variation tool
  - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.2+
  - [Bedtools](https://bedtools.readthedocs.io/en/latest/) v2
  - [Assemblytics](http://assemblytics.com/) or
  - [RaGOO](https://github.com/malonge/RaGOO)
- For Populational variation tool
  - R 3.3+ libs
    - [ggplot2](https://ggplot2.tidyverse.org/)
    - [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2=)
    - [extrafont](https://cran.r-project.org/web/packages/extrafont/README.html)
    - rmarkdown
    - kableExtra
    - dplyr
    - reshape2
    - forcats
    - ggthemes
    - rjson
    - viridisLite
    - viridis
    - bookdown
    - knitr
  - [Snakemake](https://snakemake-wrappers.readthedocs.io/en/stable/) 5.5.2+
  - [Minimap2](https://github.com/lh3/minimap2) 2.16+
  - [Samtools](http://www.htslib.org/) 1.10+
  - [Sniffles](https://github.com/fritzsedlazeck/Sniffles) 1.0.10+
  - [Flye 2.8+ - optional](https://github.com/fenderglass/Flye)
  - [WTDGB2 - optional](https://github.com/ruanjue/wtdbg2)
  - Python libs
    - [Biopython](https://biopython.org/)
    - [Pandas](https://pandas.pydata.org/)
    - [Numpy](https://numpy.org/)
    - [pylab](https://matplotlib.org/)
  - Perl v5.26.2+

## Installation<a name="Installation"></a>

### Using Git<a name="git"></a>

Once the requirements fullfilled, just git clone

```
  git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git
```

### Using Singularity<a name="singularity"></a>

#### Compiling yourself
A [*Singularity* container](https://sylabs.io/) is available with all tools compiled in. The *def* file provided can be compiled as such:

```
singularity build Singularity Singularity.TrEMOLO-2.0.def

```

**YOU MUST BE ROOT for compiling**

#### Pulling from SingularityHub

This option is disabled since Singularity Hub is for the moment in read-only


<div id='config_file'/>

# Configuration of the parameter file<a name="configuration"></a>

You will first have to enter your parameters in a *.yaml* file (see example *config.yaml* file). The necessary parameters are :

```
# all path can be relatif or absolute
DATA:
    REFERENCE:       "/path/to/reference_file.fasta"   #reference genome (fasta file) only if INSIDER_VARIANT = True [optional]
    GENOME:          "/path/to/genome_file.fasta"      #genome (fasta file) [required]
    SAMPLE:          "/path/to/reads_file.fastq"       #long reads (a fastq file) only if OUTSIDER_VARIANT = True [optional]
    WORK_DIRECTORY:  "/path/to/directory_work"         #name of output directory [required or empty]
    TE_DB:           "/path/to/database_TE.fasta"      #Database of TE (a fasta file) [required]



CHOICE:
    PIPELINE:
        OUTSIDER_VARIANT: True  # TE no assembled (out of genome)
        INSIDER_VARIANT: True   # TE assembled (in genome)
        REPORT: True            # for getting report.html with graph
    OUTSIDER_VARIANT:
        LOCAL_ASSEMBLY:
            FLYE: False # (True, False)
            WTDGB: False # (True, False)
        CALL_SV: "svim" # possibility (sniffles, svim)
        INTEGRATE_TE_TO_GENOME: True # (True, False) xxx
        OPTIMIZE_FREQUENCE: True # (True, False) xxx
    INTERMEDIATE_FILE: True     # to keep the intermediate analysis files to process them.


PARAMS:
    OUTSIDER_VARIANT:
        MINIMAP2:
            PRESET_OPTION: 'map-ont' # minimap2 preset option is map-ont by default (map-pb, map-ont etc)
            OPTION: '-t 8' # more option of minimap2
        SAMTOOLS_VIEW:
            PRESET_OPTION: ''
        SAMTOOLS_SORT:
            PRESET_OPTION: ''
        SAMTOOLS_CALLMD:
            PRESET_OPTION: ''
        FLYE:
            OPTIONS: '--plasmids -t 8' #
            PRESET_OPTION: '--nano-raw' #
        TSD:
            FILE_SIZE_TE_TSD: "/path/to/SIZE_TSD.txt"
            SIZE_FLANK: 30  # flanking sequence size to calculate TSD put value > 4
        TE_DETECTION:
            CHROM_KEEP: "." # regular expresion of chromosome; exemple  for Drosophila  "2L,2R,3[RL],X" ; Put "." for keep all chromosome
        INTEGRATE_TE_TO_GENOME:
            PUT_ID: True # (True, False) xxx
            PUT_SEQUENCE_DB_TE: True # (True, False) xxx
        PARS_BLN_OPTION: "" # option of TrEMOLO/lib/python/parse_blast_main.py d'ont put -c option
    INSIDER_VARIANT:
        PARS_BLN_OPTION: "--min-size-percent 80 --min-pident 80"


```

Main parameters

`REFERENCE` : Reference genome (fasta file) only if INSIDER_VARIANT = True
`WORK_DIRECTORY` : Directory that will contain the output files (if the directory does not exist it will be created). default value is "output".
`GENOME` : Genome assembly (fasta file)
`SAMPLE` : File containing the reads of genome assembly.
`TE_DB`  : **fasta** file containing the sequence of transposable elements.

<div id='start_pipeline'/>

# Usage<a name="usage"></a>

```
snakemake --snakefile /path/to/TrEMOLO/creation_snakefile.snk --configfile /path/to/your_config.yaml
```

after to dowload data test

```
snakemake --snakefile TrEMOLO/creation_snakefile.snk --configfile data_test_TrEMOLO/config.yml
```

# Summarize output files :open_file_folder:<a name="output"></a>

Example of output file obtained after using the pipeline.

```
WORK_DIRECTORY
├── params.log
├── params.yaml
├── POSITION_START_TE_INSIDER.bed
├── POSITION_START_TE_OUTSIDER.bed
├── FREQ_AFTER
│   ├── FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE_COUNT.csv
│   ├── FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv
│   ├── MAPPING_POSTION_TE_MD.sorted.bam
│   ├── MAPPING_POSTION_TE_MD.sorted.bam.bai
│   ├── TE_POSITION_SIZE.txt
│   ├── variants.bln
│   ├── variants.fasta
│   └── variants.vcf
├── INSIDER
│   ├── FREQ_GLOBAL
│   │   └── DEPTH_TE_INSIDER.csv
│   ├── TE_DETECTION
│   │   ├── DELETION.bed
│   │   ├── DELETION.bln
│   │   ├── DELETION_COMBINE_TE.csv
│   │   ├── DELETION_COUNT_TE.csv
│   │   ├── DELETION.csv
│   │   ├── DELETION_SEQ.fasta
│   │   ├── INSERTION.bed
│   │   ├── INSERTION.bln
│   │   ├── INSERTION_COMBINE_TE.csv
│   │   ├── INSERTION_COUNT_TE.csv
│   │   ├── INSERTION.csv
│   │   ├── INSERTION_SEQ.fasta
│   │   ├── INSERTION_TE.bed
│   │   └── TSD
│   │       ├── TSD_TE_NAME.txt
│   │       ├── TSD_17.6.txt
│   │       ├── TSD_blood.txt
│   │       ├── TSD_Idefix.txt
│   │       ├── TSD_roo.txt
│   │       ├── TSD_rover.txt
│   │       └── TSD_ZAM.txt
│   ├── TE_INSIDER_VR
│   │   ├── DELETION.bed
│   │   ├── DELETION.bln
│   │   ├── DELETION_COUNT_TE.csv
│   │   ├── DELETION.csv
│   │   ├── DELETION_SEQ.fasta
│   │   ├── INSERTION.bed
│   │   ├── INSERTION.bln
│   │   ├── INSERTION_COUNT_TE.csv
│   │   ├── INSERTION.csv
│   │   ├── INSERTION_SEQ.fasta
│   └── VARIANT_CALLING
│       ├── assemblytics_out.Assemblytics_assembly_stats.txt
│       ├── assemblytics_out.Assemblytics_structural_variants.bed
│       ├── assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta
│       ├── assemblytics_out.coords.csv
│       ├── assemblytics_out.coords.tab
│       ├── assemblytics_out.variants_between_alignments.bed
│       ├── assemblytics_out.variants_within_alignments.bed
│       ├── pm_against_ref.sam
│       └── pm_against_ref.sam.delta
├── log  #log file to check if you have any error
│   ├── DETECTION_TE.err
│   ├── extract_read.err
│   ├── extract_read.out
│   ├── FIND_TE_ON_REF.err
│   ├── FIND_TE_ON_REF.out
│   ├── FREQ_GLOBAL
│   ├── FREQUENCE.err
│   ├── FREQUENCE.out
│   ├── minimap2
│   ├── minimap2.err
│   ├── pm_contigs_against_ref.sam.log
│   ├── samtools.err
│   ├── samtools.out
│   ├── Snakefile_insider.log
│   ├── Snakefile_outsider.log
│   ├── sniffles.err
│   ├── sniffles.out
│   ├── TE_INSIDER.err
│   ├── TE_INSIDER.out
│   ├── TSD.err
│   └── TSD.out
├── OUTSIDER
│   ├── ET_FIND_FA
│   │   ├── TE_REPORT_FOUND_17.6.fasta
│   │   ├── TE_REPORT_FOUND_297.fasta
│   │   ├── TE_REPORT_FOUND_3S18.fasta
│   │   ├── TE_REPORT_FOUND_412.fasta
│   │   ├── TE_REPORT_FOUND_blood.fasta
│   │   ├── TE_REPORT_FOUND_Burdock.fasta
│   │   ├── TE_REPORT_FOUND_copia.fasta
│   │   ├── TE_REPORT_FOUND_gtwin.fasta
│   │   └── TE_REPORT_FOUND_ZAM.fasta
...
│   ├── FASTA_FIND
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta

│   ├── FIND_TE_ON_REF
│   ├── FREQ_AFTER
│   │   └── DEPTH_TE.csv
│   ├── ID_BEST_READ_TE.txt
│   ├── ID_READS_TE
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.txt
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.txt
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.txt
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.txt
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.txt

│   ├── INSIDER_VR
│   │   ├── assemblytics_out.Assemblytics_assembly_stats.txt
│   │   ├── assemblytics_out.Assemblytics_structural_variants.bed
│   │   ├── assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta
│   │   ├── assemblytics_out.coords.csv
│   │   ├── assemblytics_out.coords.tab
│   │   ├── assemblytics_out.variants_between_alignments.bed
│   │   ├── assemblytics_out.variants_within_alignments.bed
│   │   ├── pm_against_ref.sam
│   │   └── pm_against_ref.sam.delta
│   ├── MAPPING
│   │   ├── MAPPING_POSTION_TE.bam
│   │   ├── MAPPING_POSTION_TE.bam.bai
│   │   ├── SAMPLE_mapping_GENOME_MD.sorted.bam
│   │   └── SAMPLE_mapping_GENOME_MD.sorted.bam.bai
│   ├── MAPPING_TO_REF
│   │   └── INSERTION_TE.bed
│   ├── READ_FASTQ_TE
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fastq
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fastq
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fastq
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fastq
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fastq

│   ├── TE_DETECTION
│   │   ├── BLAST_SEQUENCE_INDEL_vs_DBTE.bln
│   │   ├── COMBINE_TE.csv
│   │   ├── DEPTH_TE.csv
│   │   ├── FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE_COUNT.csv
│   │   ├── FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv
│   │   ├── POSITION_START_TE.bed
│   │   └── TSD
│   │       ├── TSD_17.6_KO.txt
│   │       ├── TSD_17.6_OK.txt
│   │       ├── TSD_17.6_TSD_OK.txt
│   │       ├── TSD_17.6_TSM_OK.txt
│   │       ├── TSD_17.6.txt
│   │       ├── TSD_297.txt
│   │       ├── TSD_3S18.txt
│   │       ├── TSD_412_KO_corrected.txt
│   │       ├── TSD_412_KO.txt
│   │       ├── TSD_412_OK.txt
│   │       ├── TSD_412_TSD_OK.txt
│   │       ├── TSD_412_TSM_OK.txt
│   │       ├── TSD_412.txt
│   │       ├── TSD_blood_KO_corrected.txt
│   │       ├── TSD_blood_KO.txt
│   │       ├── TSD_blood_OK.txt
│   │       ├── TSD_blood_TSD_OK.txt
│   │       ├── TSD_blood_TSM_OK.txt
...
│   ├── TE_TOWARD_GENOME
│   │   ├── genome.out.fasta       #pseudo genome
│   │   ├── genome.out.fasta.fai
│   │   ├── TRUE_POSITION_TE.bed   #nouvelle emplacement des TE intégré
│   │   └── TRUE_POSITION_TE.fasta #sequence des TE intégré
│   └── VARIANT_CALLING
│       ├── SEQUENCE_INDEL.fasta
│       ├── SEQUENCE_INDEL.fasta.fai
│       └── SV.vcf
├── params.log
├── params.yaml
├── POSITION_TE_INSIDER.bed
├── POSITION_TE_OUTSIDER.bed
├── REPORT
│   ├── mini_report
│   └── report.html
├── SNAKE_USED
│   ├── Snakefile_insider.snk
│   └── Snakefile_outsider.snk
├── VALUES_TSD_ALL_GROUP.csv
├── VALUES_TSD_GROUP_OUTSIDER.csv
└── VALUES_TSD_INSIDER_GROUP.csv
```

The most useful output files are :

* The html report in **REPORT/report.html** with the graphs
* The positions of **TE POSITION_TE_INSIDER.bed**, **POSITION_TE_OUTSIDER.bed**

Some csv file (INSERTION.csv, FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv) has precise information on the identification of TE :

| sseqid | qseqid | pident | size_per | size_el | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore |
| ------ | ------ | ------ | -------- | ------- | -------- | ------- | ------ | ---- | ------ | ---- | ------ | -------- |
| ZAM | 2R:\<INS\>:12136769:12145149:33748:4:IMPRECISE:- | 95.494 | 99.0 | 8347 | 123 | 178 | 5 | 8352 | 8435 | 1 | 0.0 | 13369.0 |
| blood | 3R:\<INS\>:22519173:22526514:100924:1:PRECISE:+ | 94.259 | 99.0 | 7338 | 164 | 189 | 3 | 7341 | 7410 |  1 | 0.0 | 11230.0 |

## Description of header .csv file :

 1.    `qseqid` : query (e.g., gene) sequence id
 2.    `sseqid` : subject (e.g., reference genome) sequence id
 3.    `pident` : percentage of identical matches
 4.    `size_per` :   percentage of size TE
 5.    `size_el` :  size sequence TE aligned
 6.    `mismatch` : number of mismatches
 7.    `gapopen` :  number of gap openings
 8.    `qstart` :  start of alignment in query
 9.    `qend`  : end of alignment in query
 10.   `sstart`   : start of alignment in subject
 11.   `send`   : end of alignment in subject
 12.   `evalue`   : expect value
 13.   `bitscore`   : bit score

You can find the description here : [http://www.metagenomics.wiki/tools/blast/blastn-output-format-6](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)


# Licence and Citation<a name="citation"></a>

It is licencied under [CeCill-C](Licence_CeCILL-C_V1-en.txt) and [GPLv3](LICENSE).

If you use TrEMOLO, please cite:

[Mohamed, M.; Dang, N. .-M.; Ogyama, Y.; Burlet, N.; Mugat, B.; Boulesteix, M.; Mérel, V.; Veber, P.; Salces-Ortiz, J.; Severac, D.; Pélisson, A.; Vieira, C.; Sabot, F.; Fablet, M.; Chambeyron, S. A Transposon Story: From TE Content to TE Dynamic Invasion of Drosophila Genomes Using the Single-Molecule Sequencing Technology from Oxford Nanopore. Cells 2020, 9, 1776.](https://www.mdpi.com/2073-4409/9/8/1776)
