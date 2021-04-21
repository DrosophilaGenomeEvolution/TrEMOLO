
<img src="images/TrEMOLO4.png">

# TrEMOLO

*******
Content 
 1. [Configuration of the parameter file.](#config_file)
 2. [Start the pipeline.](#start_pipeline)

*******

## Require 

  * [Snakemake](https://snakemake-wrappers.readthedocs.io/en/stable/) 5.5.2
  * [Minimap2](https://github.com/lh3/minimap2) 2.16+
  * [Samtools](http://www.htslib.org/) 1.10+
  * [Sniffles](https://github.com/fritzsedlazeck/Sniffles) 1.0.10 or 1.0.12
  * svim 1.4.2
  * [Biopython](https://biopython.org/)
  * [Pandas](https://pandas.pydata.org/)
  * [Numpy](https://numpy.org/)

<div id='config_file'/> 

## Configuration of the parameter file 

You will first have to enter your parameters in a **.yaml** file (see example config.yaml file). The necessary parameters are :

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
            FLYE: False
            WTDGB: False
        CALL_SV: "svim" # possibility (sniffles, svim) 
        INTEGRATE_TE_TO_GENOME: True
        OPTIMIZE_FREQUENCE: True
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
            CHROM_KEEP: "." # regular expresion of chromosome; exemple  for Drosophila  "2L,2R,3[RL],X" ; put "." for keep all chromosome
        INTEGRATE_TE_TO_GENOME:
            PUT_ID: True
            PUT_SEQUENCE_DB_TE: True
        PARS_BLN_OPTION: "" # option of TrEMOLO/pipeline/lib/python/parse_blast_main.py d'ont put -c option
    INSIDER_VARIANT:
        PARS_BLN_OPTION: "--min-size-percent 80 --min-pident 80" 
        
        
```

`WORK_DIRECTORY` : directory that will contain the output files (if the directory does not exist it will be created).
`GENOME` : Genome assembly
`SAMPLE` : File containing the reads of genome assembly.
`TE_DB` : File **.fasta** containing the sequence of transposable elements.

<div id='start_pipeline'/> 

## Start the pipeline

```
snakemake --snakefile /path/to/TrEMOLO/pipeline/creation_snakefile.snk --configfile /path/to/your_config.yaml
```


### Summarize output files :open_file_folder:

Example of output file obtained after using the pipeline, in work directory.

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
│   ├── tmp_MAPPING_POSTION_TE_MD.sorted.bam
│   ├── tmp_MAPPING_POSTION_TE.sam
│   ├── tmp_MAPPING_POSTION_TE.sorted.bam
│   ├── tmp_MAPPING_POSTION_TE.sorted.bam.bai
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
│   │       ├── TSD_1360.txt
│   │       ├── TSD_17.6.txt
│   │       ├── TSD_297.txt
│   │       ├── TSD_3S18.txt
│   │       ├── TSD_412.txt
│   │       ├── TSD_Bari1.txt
│   │       ├── TSD_blood.txt
│   │       ├── TSD_BS.txt
│   │       ├── TSD_Burdock.txt
│   │       ├── TSD_copia.txt
│   │       ├── TSD_diver.txt
│   │       ├── TSD_Doc.txt
│   │       ├── TSD_FB.txt
│   │       ├── TSD_F-element.txt
│   │       ├── TSD_flea.txt
│   │       ├── TSD_G2.txt
│   │       ├── TSD_gypsy5.txt
│   │       ├── TSD_gypsy6.txt
│   │       ├── TSD_HB.txt
│   │       ├── TSD_HMS-Beagle2.txt
│   │       ├── TSD_HMS-Beagle.txt
│   │       ├── TSD_hobo.txt
│   │       ├── TSD_hopper.txt
│   │       ├── TSD_Idefix.txt
│   │       ├── TSD_I-element.txt
│   │       ├── TSD_invader3.txt
│   │       ├── TSD_Ivk.txt
│   │       ├── TSD_jockey.txt
│   │       ├── TSD_Juan.txt
│   │       ├── TSD_Max-element.txt
│   │       ├── TSD_McClintock.txt
│   │       ├── TSD_mdg1.txt
│   │       ├── TSD_mdg3.txt
│   │       ├── TSD_micropia.txt
│   │       ├── TSD_opus.txt
│   │       ├── TSD_P-element.txt
│   │       ├── TSD_pogo.txt
│   │       ├── TSD_Quasimodo.txt
│   │       ├── TSD_roo.txt
│   │       ├── TSD_rover.txt
│   │       ├── TSD_Rt1a.txt
│   │       ├── TSD_Rt1b.txt
│   │       ├── TSD_S-element.txt
│   │       ├── TSD_springer.txt
│   │       ├── TSD_Stalker2.txt
│   │       ├── TSD_Stalker4.txt
│   │       ├── TSD_Tc1.txt
│   │       ├── TSD_Transpac.txt
│   │       ├── TSD_X-element.txt
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
│   │   └── tmp.fasta
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
├── log
│   ├── DETECTION_TE
│   ├── DETECTION_TE.err
│   ├── extract_read
│   ├── extract_read.err
│   ├── extract_read.out
│   ├── FIND_TE_ON_REF
│   ├── FIND_TE_ON_REF.err
│   ├── FIND_TE_ON_REF.out
│   ├── FREQ_GLOBAL
│   ├── FREQUENCE
│   ├── FREQUENCE.err
│   ├── FREQUENCE.out
│   ├── minimap2
│   ├── minimap2.err
│   ├── pm_contigs_against_ref.sam.log
│   ├── report
│   ├── samtools
│   ├── samtools.err
│   ├── samtools.out
│   ├── Snakefile_insider.log
│   ├── Snakefile_outsider.log
│   ├── sniffles
│   ├── sniffles.err
│   ├── sniffles.out
│   ├── TE_INSIDER
│   ├── TE_INSIDER.err
│   ├── TE_INSIDER.out
│   ├── TE_TOWARD_GENOME
│   ├── TSD
│   ├── TSD.err
│   └── TSD.out
├── OUTSIDER
│   ├── ET_FIND_FA
│   │   ├── TE_REPORT_find_17.6.fasta
│   │   ├── TE_REPORT_find_297.fasta
│   │   ├── TE_REPORT_find_3S18.fasta
│   │   ├── TE_REPORT_find_412.fasta
│   │   ├── TE_REPORT_find_blood.fasta
│   │   ├── TE_REPORT_find_Burdock.fasta
│   │   ├── TE_REPORT_find_copia.fasta
│   │   ├── TE_REPORT_find_diver.fasta
│   │   ├── TE_REPORT_find_Doc.fasta
│   │   ├── TE_REPORT_find_F-element.fasta
│   │   ├── TE_REPORT_find_flea.fasta
│   │   ├── TE_REPORT_find_G2.fasta
│   │   ├── TE_REPORT_find_gtwin.fasta
│   │   ├── TE_REPORT_find_gypsy6.fasta
│   │   ├── TE_REPORT_find_gypsy.fasta
│   │   ├── TE_REPORT_find_HMS-Beagle.fasta
│   │   ├── TE_REPORT_find_Max-element.fasta
│   │   ├── TE_REPORT_find_mdg1.fasta
│   │   ├── TE_REPORT_find_mdg3.fasta
│   │   ├── TE_REPORT_find_NOF.fasta
│   │   ├── TE_REPORT_find_opus.fasta
│   │   ├── TE_REPORT_find_Quasimodo.fasta
│   │   ├── TE_REPORT_find_roo.fasta
│   │   ├── TE_REPORT_find_rover.fasta
│   │   ├── TE_REPORT_find_springer.fasta
│   │   ├── TE_REPORT_find_Stalker4.fasta
│   │   ├── TE_REPORT_find_Transpac.fasta
│   │   └── TE_REPORT_find_ZAM.fasta
│   ├── FASTA_FIND
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta.nin
│   │   ├── reads_2L_RaGOO_RaGOO:100:1621121-1621132.fasta.nsq
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta.nin
│   │   ├── reads_2L_RaGOO_RaGOO:122:1727638-1727638.fasta.nsq
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fasta.nin
│   │   ├── reads_2L_RaGOO_RaGOO:145:2570185-2570185.fasta.nsq
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fasta.nin
│   │   ├── reads_2L_RaGOO_RaGOO:166:4083207-4083207.fasta.nsq
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fasta.nin
│   │   ├── reads_2L_RaGOO_RaGOO:167:4091408-4091408.fasta.nsq
│   │   ├── reads_2L_RaGOO_RaGOO:168:4096613-4096613.fasta
│   │   ├── reads_2L_RaGOO_RaGOO:168:4096613-4096613.fasta.fai
│   │   ├── reads_2L_RaGOO_RaGOO:168:4096613-4096613.fasta.nhr
│   │   ├── reads_2L_RaGOO_RaGOO:168:4096613-4096613.fasta.nin

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
│   │   ├── FLANK_TE.bed
│   │   ├── FLANK_TE_PASS2.bed
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
│   │       ├── DIR_SEQ_TE_READ_POS
│   │       │   ├── sequence_TE_reads_2L_RaGOO_RaGOO:100:1621121-1621132..bed
│   │       │   ├── sequence_TE_reads_2L_RaGOO_RaGOO:100:1621121-1621132..bln
│   │       │   ├── sequence_TE_reads_2L_RaGOO_RaGOO:122:1727638-1727638..bed
│   │       │   ├── sequence_TE_reads_2L_RaGOO_RaGOO:122:1727638-1727638..bln
│   │       │   ├── sequence_TE_reads_2L_RaGOO_RaGOO:145:2570185-2570185..bed

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

│   ├── TE_TOWARD_GENOME
│   │   ├── genome.out.fasta
│   │   ├── genome.out.fasta.fai
│   │   ├── ID.txt
│   │   ├── SEQUENCE_INDEL_as_TE.bed
│   │   ├── SEQUENCE_INDEL_as_TE.fasta
│   │   ├── SEQUENCE_INDEL.bed
│   │   ├── SEQUENCE_INDEL_DB_TE.fasta
│   │   ├── SEQUENCE_INDEL_TE.bed
│   │   ├── TRUE_POSITION_TE.bed
│   │   └── TRUE_POSITION_TE.fasta
│   └── VARIANT_CALLING
│       ├── SEQUENCE_INDEL.fasta
│       ├── SEQUENCE_INDEL.fasta.fai
│       └── SV.vcf
├── params.log
├── params.yaml
├── POSITION_START_TE_INSIDER.bed
├── POSITION_START_TE_OUTSIDER.bed
├── REPORT
│   ├── mini_report
│   └── report.html
├── SNAKE_USED
│   ├── Snakefile_insider.snk
│   └── Snakefile_outsider.snk
├── VALUES_TSD_ALL_GROUP.csv
├── VALUES_TSD_GROUP.csv
├── VALUES_TSD_INSIDER.csv
└── VALUES_TSD_INSIDER_GROUP.csv
```

In directory **ET_FIND_FA** you can found fasta files of all sequences TE by family of TE


The output file whose name ends with "cnTE_ALL_ET.csv" contains the following informations :

| sseqid | qseqid | pident | size_per | size_el | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore |
| ------ | ------ | ------ | -------- | ------- | -------- | ------- | ------ | ---- | ------ | ---- | ------ | -------- |
| ZAM | 2R:\<INS\>:12136769:12145149:33748:4:IMPRECISE | 95.494 | 99.0 | 8347 | 123 | 178 | 5 | 8352 | 8435 | 1 | 0.0 | 13369.0 |
| blood | 3R:\<INS\>:22519173:22526514:100924:1:PRECISE | 94.259 | 99.0 | 7338 | 164 | 189 | 3 | 7341 | 7410 |  1 | 0.0 | 11230.0 |

##### Description of header .csv file :

 1.    `qseqid` :   query (e.g., gene) sequence id
 2.    `sseqid` : subject (e.g., reference genome) sequence id
 3.    `pident` :  percentage of identical matches
 4.    `size_per`    :   percentage of size TE
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

##### Description of format qseqid or header of fasta files  :

chromosome:type_variant:start_position:end_position:ID:NB_read_support:PRECISE_OR_IMPRECISE_POSITION_SV
3R:\<INS\>:22519173:22526514:100924:1:PRECISE

 1.    `chromosome` : the name of the chromosome
 2.    `type_variant` : type of variant (insertion, deletion, duplication...) provided by sniffle/svim in the **vcf** file
 3.    `start_position` : 5' position of the TE on the chromosome
 4.    `end_position`  : 3' position of the TE on the chromosome
 5.    `ID` : a uniq id of variant provided by sniffle in the **vcf** file
 6.    `NB_read_support` : number of reads supports
 7.    `PRECISE_OR_IMPRECISE_POSITION_SV` : indicates if the location of the TE is precise or imprecise


Note : if you want to use the scripts individually know that the format **qseqid** is a format necessary for most scripts in the pipeline. The header of some fasta files must also be in this format. This info is extracted from a vcf file (generate by sniffles/svim with the appropriate options).

In directory **FASTA_FIND**, **REGION_RD_G73vsG73LR**, **READ_FASTQ_G73vsG73LR** the format of name file is : 

reads_chromosome:ID:start_position:end_position.(txt|fastq|fasta)


##### Diagram of the different stages (rules)

<img src="images/dag.svg">

Enjoy :+1:


