# TrEMOLO
Transposable Elements MOvement detection using LOng reads

TrEMOLO used long reads assemblies to detect
- Global TE variations between two assembled genomes  
- Populational variation in TE insertions.

It is licencied under [CeCill-C](Licence_CeCILL-C_V1-en.txt) and [GPLv3](LICENSE).

Please cite xxxx if using TrEMOLO

### Requirements:

- For both approaches
  - Python 3.6+
- For Global variation tool
  - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.2+
  - [Bedtools](https://bedtools.readthedocs.io/en/latest/) v2
  - [Assemblytics](http://assemblytics.com/) or [RaGOO](https://github.com/malonge/RaGOO) output
- For Populational variation tool
  - [Snakemake](https://snakemake-wrappers.readthedocs.io/en/stable/) 5.5.2
  - [Minimap2](https://github.com/lh3/minimap2) 2.16+
  - [Samtools](http://www.htslib.org/) 1.10+
  - [Sniffles](https://github.com/fritzsedlazeck/Sniffles) 1.0.10+
  - Biopython
  - Pandas
  - Numpy

### Installation:
Once the requirements fullfilled, just git clone

```
  git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git
  cd TrEMOLO
```

## Global TE variations
The *svTEidentification.py* script will inform about the new TE putative insertion based on a Assemblytics BED file. This file ca be created through Assemblytics or RaGOO.
  A TE database (mutlifasta file) must be provided. You can format it using the *makeblastdb* command from BLAST suite.

  Requires Bedtools and BLAST to be accessible in the path

```
 svTEidentification.py -i Assemblytics.bed -d TEdatabase -o output
```
Options
```
  - \-h, --help
          show this help message and exit
  - \-v, --version
          display svTEidentification.py version number and exit

Input mandatory infos for running:
  - \-i <filename>, --input <filename>
          BED file issued from Assemblytics/RaGOO analysis
  - \-r <filename>, --reference <filename>
          reference sequence file used in Assemblytics/RaGOO
      - \-a <filename>, --alternate <filename>
          alternate sequence file used in Assemblytics/RaGOO
      - \- d <filename>, --database <filename>
          TE database in multifasta format (must be formatted)
  - \-o <filename>, --out <filename>
          Prefix of output files
  - \-s <minimalPercentage>, --size <minimalPercentage>
          Minimal percentage of identity and size for a TE hit to be conserved (Optional, default 90)
```
This script will output 8 files:
- Two BED files representing respectively the deletions and insertions of the Alternate relative to the Reference
- Two FASTA files with the sequences corresponding to these deletions and insertions
- Two BLASTn output files providing the whole results of the precedent files compared to the TE database
- And finally, two tabulated files (.csv) coming from filtering the BLASTn results according to the minimal percentage threshold.

The tabulated files are structured as follows:


|#TE | Location | PercId | FragSize | RefSize | PercTotal |
| ------------- | ------------- | ------------- | ------------- |------------- |------------- |
| blood | +::2R_RaGOO_RaGOO:9104794-9112210 | 99.569 | 7424 | 7410 | 100.2 |
| 412 | +::2R_RaGOO_RaGOO:13514820-13522246 | 99.511 | 6134 | 7567 | 81.06 |
| flea | +::X_RaGOO_RaGOO:7469447-7474489 | 99.543 | 5036 | 5034 | 100.04 |
