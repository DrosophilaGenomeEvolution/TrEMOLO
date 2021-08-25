

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5391)


<img src="images/TrEMOLO9.png">

- [Introduction](#introduction)
    - [Global variations](#in)
    - [Populational variations](#out)
- [Release note](#release)
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

## Global variations, the insiders<a name="in"></a>

Using a reference genome and an assembled one (preferentially using long contigs or even better a chrosomome-scale assembly), TrEMOLO will extract the **insiders**, *i.e.* variant transposable elements (TEs) present globally in the assembly, and tag them. Indeed, assemblers will provide the most frequent haplotype at each locus, and thus an assembly represent just the "consensus" of all haplotypes present at each locus.
You will obtain a [set of files](#output) with the location of these variable insertions.

## Populational variations, the outsiders<a name="out"></a>

Through remapping of reads that have been used to assemble the genome of interest, TrEMOLO will identify the populational variations (or even somatic ones) within the initial dataset of reads, and thus of DNA/individuals sampled. These variant TEs are the **outsiders**, present only in a subpart of the population.
In the same way as for insiders, you will obtain a [set of files](#output) with the location of these variable insertions.


# Release notes<a name="release"></a>

# Requirements<a name="requirements"></a>

Numerous tools are used by TrEMOLO. We recommand to use the [Singularity installation](#singularity) to be sure to have all of them in the good configuration.

- For both approaches
  - Python 3.6+
- For Global variation tool
  - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.2+
  - [Bedtools 2.27.1](https://bedtools.readthedocs.io/en/latest/) v2
  - [Assemblytics](http://assemblytics.com/) or
  - [RaGOO](https://github.com/malonge/RaGOO)
- For Populational variation tool
  - [Snakemake](https://snakemake-wrappers.readthedocs.io/en/stable/) 5.5.2+
  - [Minimap2](https://github.com/lh3/minimap2) 2.16+
  - [Samtools](http://www.htslib.org/) 1.9+
  - [Sniffles 1.0.12](https://github.com/fritzsedlazeck/Sniffles)
  - [Flye 2.8+ - optional](https://github.com/fenderglass/Flye)
  - [WTDGB2 - optional](https://github.com/ruanjue/wtdbg2)
  - pandoc-citeproc 0.17
  - Python libs
    - [Biopython](https://biopython.org/)
    - [Pandas](https://pandas.pydata.org/)
    - [Numpy](https://numpy.org/)
    - [pylab](https://matplotlib.org/)
    - [intervaltree](https://pypi.org/project/intervaltree/)
  - Perl v5.26.2+
- For report
  - R 3.3+ libs
    - [ggplot2](https://ggplot2.tidyverse.org/)
    - [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2=)
    - [extrafont](https://cran.r-project.org/web/packages/extrafont/README.html)
    - [rmarkdown](https://rmarkdown.rstudio.com/)
    - [kableExtra](https://bookdown.org/yihui/rmarkdown-cookbook/kableextra.html)
    - [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
    - [reshape2](https://www.r-project.org/nosvn/pandoc/dplyr.html)
    - [forcats](https://rdrr.io/cran/forcats/)
    - [ggthemes](https://github.com/jrnold/ggthemes)
    - [rjson](https://rdrr.io/cran/rjson/)
    - [viridisLite](https://github.com/sjmgarnier/viridisLite)
    - [viridis](https://www.rdocumentation.org/packages/viridis/versions/0.3.4)
    - [bookdown](https://bookdown.org/yihui/bookdown/get-started.html)
    - [knitr](https://www.r-project.org/nosvn/pandoc/knitr.html)

# Installation<a name="Installation"></a>

## Using Git<a name="git"></a>

Once the requirements fullfilled, just git clone

```bash
git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git
```

## Using Singularity<a name="singularity"></a>


[*Singularity* installation Debian/Ubuntu with package](https://sylabs.io/guides/3.0/user-guide/installation.html#install-the-debian-ubuntu-package-using-apt)

### Compiling yourself
A [*Singularity* container](https://sylabs.io/) is available with all tools compiled in.
The *def* file provided can be compiled as such:

```bash
sudo singularity build TrEMOLO.simg TrEMOLO/Singularity
```

Test TrEMOLO with singularity

```bash
singularity exec TrEMOLO.simg snakemake --snakefile TrEMOLO/creation_snakefile.snk --configfile TrEMOLO/test/tmp_config.yml
```

**YOU MUST BE ROOT for compiling**

### Pulling from SingularityHub

This option is disabled since Singularity Hub is for the moment in read-only



# Configuration of the parameter file<a name="configuration"></a>

TrEMOLO uses snakemake to perform its analyses. You have then first to provide your parameters in a *.yaml* file (see an example in the *config.yaml* file). Parameters are :

```yaml
# all path can be relative or absolute depending of your tree. It is advised to only use absolute path if you are not familiar with computer science or the importance of folder trees structure. 
DATA:
    GENOME:          "/path/to/genome_file.fasta"      #genome (fasta file) [required]
    TE_DB:           "/path/to/database_TE.fasta"      #Database of TE (a fasta file) [required]
    REFERENCE:       "/path/to/reference_file.fasta"   #reference genome (fasta file) only if INSIDER_VARIANT = True [optional]
    SAMPLE:          "/path/to/reads_file.fastq"       #long reads (a fastq file) only if OUTSIDER_VARIANT = True [optional]
    #At least, provide either REFERENCE or SAMPLE. Both can be provided
    WORK_DIRECTORY:  "/path/to/directory"         #name of output directory [optional, will be created as 'TrEMOLO_OUTPUT']


CHOICE:
    PIPELINE:
        OUTSIDER_VARIANT: True  # outsiders, TE not assembled (out of assembly)
        INSIDER_VARIANT: True   # insiders, TE assembled (in assembly)
        REPORT: True            # for getting a report.html file with graphics
    OUTSIDER_VARIANT:
        CALL_SV: "svim" # possibilities: sniffles, svim
        INTEGRATE_TE_TO_GENOME: True # (True, False) Re-build the assembly with insiders integrated in
        OPTIMIZE_FREQUENCE: True # (True, False) xxx
    INTERMEDIATE_FILE: True     # to conserve the intermediate analyses files to process them.


PARAMS:
    OUTSIDER_VARIANT:
        MINIMAP2:
            PRESET_OPTION: 'map-ont' # minimap2 option is map-ont by default (map-pb, map-ont)
            OPTION: '-t 8' # more option of minimap2 can be trasnfered here
        SAMTOOLS_VIEW:
            PRESET_OPTION: ''
        SAMTOOLS_SORT:
            PRESET_OPTION: ''
        SAMTOOLS_CALLMD:
            PRESET_OPTION: ''
        TSD:
            FILE_SIZE_TE_TSD: "/path/to/SIZE_TSD.txt" # File of TSD sizes for reference elements (format="TE SIZE" by line)
            SIZE_FLANK: 30  # flanking sequence size for calculation of TSD; put value > 4
        TE_DETECTION:
            CHROM_KEEP: "." # regular expresion for chromosome filtering; for instance for Drosophila  "2L,2R,3[RL],X" ; Put "." to keep all chromosome
            GET_SEQ_REPORT_OPTION: "-m 1000" #option get_seq_vcf.py option of the sequence recovery file in the vcf
        PARS_BLN_OPTION: "--min-size-percent 90 --min-pident 94" # option for TrEMOLO/lib/python/parse_blast_main.py - don't put -c option
    INSIDER_VARIANT:
        PARS_BLN_OPTION: "--min-size-percent 80 --min-pident 80" # parameters for validation of insiders


```

The main parameters are:

- `GENOME` : Assembly of the sample of interest (or mix of samples), fasta file.\
- `TE_DB`  : **Multifasta** file containing the canonical sequence of transposable elements. You can add also copy sequences but results will be more complex to interpretate. Make sure the headers of this multifasta do not contains this caracters : ":", "|" or "/".\
- `REFERENCE` : Fasta file containing the reference genome of the sample of interest.\
- `WORK_DIRECTORY` : Directory that will contain the output files. If the directory does not exist it will be created;  default value is *output*.\
- `SAMPLE` : File containing the reads used for the sample assembly.\


You can take **config_INSIDER.yaml** for only **INSIDER** analysis or **config_OUTSIDER.yaml** for only **OUTSIDER** analysis.
To analyse **INSIDER**, only the `REFERENCE` , the `GENOME`, the `TE_DB` and the `WORK_DIRECTORY` are required. 
To analyse **OUTSIDER**, only the `SAMPLE` , the `GENOME`, the `TE_DB` and the `WORK_DIRECTORY` are required. 


# Usage<a name="usage"></a>

```bash
snakemake --snakefile /path/to/TrEMOLO/creation_snakefile.snk --configfile /path/to/your_config.yaml
```

Run test

```bash
snakemake --snakefile TrEMOLO/creation_snakefile.snk --configfile TrEMOLO/test/tmp_config.yml
```


# Output files summary :open_file_folder:<a name="output"></a>

Here is the structure of the output files obtained after running the pipeline.

```
WORK_DIRECTORY
├── params.yaml
├── POS_ALT_TE_OUTSIDER_ON_REF.bed -> INSIDER/TE_INSIDER_VR/POS_ALT_TE_OUTSIDER_ON_REF.bed
├── POSITION_TE_INOUTSIDER.bed
├── POSITION_TE_INSIDER.bed
├── POSITION_TE_OUTSIDER.bed
├── POS_TE_INSIDER_ON_REF.bed -> INSIDER/TE_DETECTION/INSERTION_TE_ON_REF.bed #POSITION TE INSIDER ON REFRENCE GENOME
├── POS_TE_OUTSIDER_ON_REF.bed #POSITION TE OUTSIDER ON REFRENCE GENOME
├── POSITION_TE_OUTSIDER_IN_NEO_GENOME.bed  #POSITION TE SEQUENCE ON BEST READS SUPPORT INTEGRATED IN GENOME
├── POSITION_TE_ON_REF.bed  #POSITION TE INSIDER AND OUTSIDER ON REFRENCE GENOME
├── POSITION_TE_OUTSIDER_IN_PSEUDO_GENOME.bed
├── VALUES_TSD_ALL_GROUP.csv
├── VALUES_TSD_GROUP_OUTSIDER.csv
├── VALUES_TSD_INSIDER_GROUP.csv
├── TE_INFO.csv
├── FREQ_OPTIMIZED
├── INSIDER #FOLDER CONTAINS FILES TRAITEMENT INSIDER
│   ├── FREQ_GLOBAL 
│   ├── TE_DETECTION
│   │   └── TSD
│   ├── TE_INSIDER_VR
│   └── VARIANT_CALLING
├── log  #log file to check if you have any error
├── OUTSIDER
│   ├── ET_FIND_FA
│   │   ├── TE_REPORT_FOUND_TE_NAME.fasta
│   │   ├── TE_REPORT_FOUND_blood.fasta
│   │   └── TE_REPORT_FOUND_ZAM.fasta
...
│   ├── FIND_TE_ON_REF
│   ├── FREQ_AFTER
│   │   └── DEPTH_TE.csv
│   ├── INSIDER_VR
│   ├── MAPPING #FOLDER CONTAINS FILES MAPPING ON GENOME
│   ├── MAPPING_TO_REF #FOLDER CONTAINS FILES MAPPING ON REFERENCE GENOME
│   ├── READ_FASTQ_TE #FOLDER CONTAINS ALL THE READs ASSOCIATED WITH THE TE
│   ├── TE_DETECTION
│   │   └── TSD
│   ├── TE_TOWARD_GENOME #FOLDER CONTAINS ALL THE READs ASSOCIATED WITH THE TE
│   │   ├── NEO_GENOME.fasta   #GENOME CONTAINS TE OUTSIDER (the best sequence of svim/sniffles)
│   │   ├── PSEUDO_GENOME_TE_DB_ID.fasta   #GENOME CONTAINS TE OUTSIDER (the sequence of database TE and the ID of svim/sniffles)
│   │   ├── TRUE_POSITION_TE.bed   #POSITION IN PSEUDO GENOME
│   │   ├── TRUE_POSITION_TE.fasta  #SEQUENCE INTEGRATE IN PSEUDO GENOME
│   │   ├── TRUE_POSITION_TE_READS.bed  #POSITION IN NEO GENOME
│   │   └── TRUE_POSITION_TE_READS.fasta  #SEQUENCE INTEGRATE IN NEO GENOME
│   └── VARIANT_CALLING  #FOLDER CONTAINS FILES OF sniflles/svim
├── REPORT
│   ├── mini_report
│   └── report.html
├── SNAKE_USED
│   ├── Snakefile_insider.snk
└── └── Snakefile_outsider.snk
```

### Most useful output

The most useful output files are :

* The html report in **your_work_direcetory/REPORT/report.html** with summary graphics, exemple [here](https://raw.githack.com/DrosophilaGenomeEvolution/TrEMOLO/v2/test/web/index.html)

The output file **your_work_direcetory/TE_INFO.csv** gathers all the necessary information.

|      chrom      |  start   | end      |   TE\|ID   |   strand  |    TSD   | SIZE_TE |      NEW_POS     |  FREQ (%) | FREQ_OPTIMIZED (%) | ID_TrEMOLO  |
| --------------- | -------- | -------- | ---------- | --------- | -------- | ------- | ---------------- | ------- | -------------- | -------------- |
|  2R_RaGOO_RaGOO | 16943971 | 16943972 | roo\|svim.INS.175 |     +     |   GTACA   | 9006    | DEFAULT:16943971 | 28.5714 |    28.5714     |    TE_ID_OUTSIDER.94047.INS.107508.0  |
|  X_RaGOO_RaGOO  | 21629415 | 21629416 | ZAM\|Assemblytics_w_534  |     -     |   CGCG   | 8435    | 21629413         | 11.1111 |    10.0000     |    INSIDER TE_ID_INSIDER.77237.Repeat_expansion.8 |


 1.    `chrom` : chromosome
 2.    `start` : start TE
 3.    `end` : end TE
 4.    `TE|ID` :   TE name and ID in **SV.vcf** (for OUTSIDER) or **assemblytics_out.Assemblytics_structural_variants.bed** (for INSIDER)
 5.    `strand` :  strand TE
 6.    `TSD` : TSD SEQUENCE 
 7.    `SIZE_TE` :  TE size
 8.    `NEW_POS` :  position corrected iwith calculate TSD (only for OUTSIDER) 
 9.    `FREQ`  : frequence normal
 10.   `FREQ_OPTIMIZED`  : frequence optimized with conversion of clipped read to not clipped
 11.   `ID_TrEMOLO`  : TrEMOLO ID


# Licence and Citation<a name="citation"></a>

It is licencied under [CeCill-C](Licence_CeCILL-C_V1-en.txt) and [GPLv3](LICENSE).

If you use TrEMOLO, please cite:

[Mohamed, M.; Dang, N. .-M.; Ogyama, Y.; Burlet, N.; Mugat, B.; Boulesteix, M.; Mérel, V.; Veber, P.; Salces-Ortiz, J.; Severac, D.; Pélisson, A.; Vieira, C.; Sabot, F.; Fablet, M.; Chambeyron, S. A Transposon Story: From TE Content to TE Dynamic Invasion of Drosophila Genomes Using the Single-Molecule Sequencing Technology from Oxford Nanopore. Cells 2020, 9, 1776.](https://www.mdpi.com/2073-4409/9/8/1776)
