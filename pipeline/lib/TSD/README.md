# TSD

*******
Content 
 1. [TSD command](#tsd)
 2. [Output](#output)

*******

## Require 
  - Python 3.6+
  - Python libs
    - [Biopython](https://biopython.org/)
    - [Pandas](https://pandas.pydata.org/)
    - [Numpy](https://numpy.org/)


<div id='tsd'/> 


## TSD command


```
usage : find_fq_to_fasta.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <NAME_OUT_DIRECTORY>
```

``find_fq_to_fasta.sh`` converted fastq files to fasta files


the script **find_fq_to_fasta.sh** requires the presence of the script **fastq_to_fasta.py** in the same directory

```
sh find_fq_to_fasta.sh work_directory/ET_FIND_FA/G73vsG73TEST_find_ZAM.fasta work_directory/READ_FASTQ_G73vsG73TEST/ FASTA_FIND ;
```

``tsd_te.sh`` search for tsd of TE


the script **tsd_te.sh** requires the presence of the script **find_tsd.py** script in the same directory

```
usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size>
```

```
sh tsd_te.sh work_directory/ET_FIND_FA/G73vsG73TEST_find_ZAM.fasta READ_FASTQ_G73vsG73TEST FASTA_FIND cannonical_TE.fasta 10 4 ;
```

<div id='output'/> 

## Output

in your results of TSD files you have for each TE 5 lines

Exemple :

```
work_G73vsG73TEST/FASTA_FIND/reads_2R_RaGOO_RaGOO:33748:12136769-12145149.fasta
>2R_RaGOO_RaGOO:<INS>:12136769:12145149:33748:4:IMPRECISE:4-8352
(CGCG, CGCG, [OK:33748], 6, 0, 4)
>CGAGTGCGCG:CGCGCTAACT:6:0:4:10:(+)
CGAGTG++:CGCG:++--|AGTTACCGACCTG...|--++:CGCG:++CTAACT
```

Description :

1. `work_directory/FASTA_FIND/reads_chrom:ID_SV:START_POSITION-END_POSITION.fasta` : the fasta file where you can find the sequence of read. 
2. `>chromosome:type_variant:start_position:end_position:ID:NB_read_support:PRECISE_OR_IMPRECISE_POSITION_SV:start_position_in_sequence_report-end_position_in_sequence_report`

	1.    `chromosome` : the name of the chromosome
	2.    `type_variant` : type of variant (insertion, deletion, duplication...) provided by sniffle
	3.    `start_position` : 5' position of the TE on the chromosome
	4.    `end_position`  : 3' position of the TE on the chromosome
	5.    `ID` : a uniq id of variant provided by sniffle
	6.    `NB_read_support` : number of reads supports
	7.    `PRECISE_OR_IMPRECISE_POSITION_SV` : indicates if the location of the TE is precise or imprecise

3. `(TSD_LEFT, TSD_RIGHT, [OK|K-O|KO:ID_SV_IN_VCF_FILE], POSITION_TSD_ON_FLANKING_SEQUENCE_LEFT, POSITION_TSD_ON_FLANKING_SEQUENCE_RIGHT, TSD_SIZE)`
	OK  : TSD located exactly at the 5' and 3' terminal of the TE and with 100% homology
	K-O : TSD located exactly at terminal 5' and 3' of the TE with 1 authorized mismatch 
	KO  : TSD not located exactly at terminal 5' and 3' but with 100% homology

4. `>FLANKING_SEQUENCE_LEFT:FLANKING_SEQUENCE_RIGHT:POSITION_TSD_ON_FLANKING_SEQUENCE_LEFT:POSITION_TSD_ON_FLANKING_SEQUENCE_RIGHT:TSD_SIZE:FLANKING_SIZE`
5. `FLANKING_SEQUENCE_LEFT++:TSD_LEFT:++--|TE_SEQUENCE|--++:TSD_RIGHT:++FLANKING_SEQUENCE_RIGHT`

A summary of the results is at the end of the file :

Exemple :

```
OK/total : 25/51
KO/total : 23/51
OK+KO/total : 48/51
K-O/total : 3/51
OK+K-O/total : 28/51
OK% : 49%
```

Enjoy :+1:


