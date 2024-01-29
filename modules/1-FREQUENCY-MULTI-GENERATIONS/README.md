# MODULE SCATTER FREQUENCY TE TrEMOLO

This graph is particularly useful for researchers aiming to track the evolution of the frequency of germinal insertions of transposable elements (TEs) across generations. It provides a clear visualization of the dynamics of these elements within the genome and offers valuable knowledge about their behavior and potential for adaptation or evolutionary change within populations over extended periods.

## Input Preparation

To function, this module requires a file containing paths to various analyses performed by TrEMOLO (minimum version +v2.5.1).

Example:

```
/path/to/work_directory_1
/path/to/work_directory_2
/path/to/work_directory_3
```

The order of paths is important as it designates the sequence in which the module will consider generations from oldest to newest.

An alternative solution allows you to manually specify the timing of generations (format work_directory:G[NUMBER]), for 

Example:

```
work_directory_1:G1
work_directory_2:G10
work_directory_3:G3
```

**Recommendation** : Do not exceed 7 paths.

This format indicates the order and the time gap between generations. In the above example, `G1` is the oldest generation, while `G10` is the most recent.

**Info** : The generation numbers (e.g., G2, G10...) enable the module to identify which Transposable Elements (TE) increase, decrease, or vary over generations.


## Run Build Graph

Execute the module with the following command:

```
singularity exec TrEMOLO.simg TrEMOLO/modules/1-FREQUENCY-MULTI-GENERATIONS/buildFrequencyGenerations.sh -i <input-init-file> [-o OUTPUT-NAME-DIRECTORY] [-g GENOME-FASTA-FILE] [-c REGEX-CHROM]
```

* `-i <input-init-file>` (required): This is the file containing the paths to your work directories.
* `[-o OUTPUT-NAME-DIRECTORY]`: Specify the name of the output directory.
* `[-g GENOME-FASTA-FILE]`: The genome (.fasta file) that was used across all work directories.
* `[-c REGEX-CHROM]`: Use this option if you want to filter by chromosome in the TE_INFOS.bed file.

For the module to function:
* Each work_directory must contain a `TE_INFOS.bed` file (output of TrEMOLO).
* If the genome is not passed with the -g option, `work_directory/OUTSIDER/FREQUENCY/MAPPING_POSITION_TE.bam` files are necessary.
* The same genome (GENOME parameter in TrEMOLO) must be used for all `work_directories`.

If you wish to select specific TEs, create a file `work_directory/TE_FREQUENCY_TrEMOLO.bed` in the same format as `TE_INFOS.bed`, including only the `OUTSIDER`.

The `[-c REGEX-CHROM]` option is useless if you have put `work_directory/TE_FREQUENCY_TrEMOLO.bed` files in your work directories.

For running tests

```bash
singularity exec TrEMOLO.simg TrEMOLO/modules/1-FREQUENCY-MULTI-GENERATIONS/buildFrequencyGenerations.sh \
    -i TrEMOLO/modules/1-FREQUENCY-MULTI-GENERATIONS/test/INIT_FREQ_TE_TrEMOLO.txt \
    -o TEST-GRAPH-FRQUENCIES \
    -c "^[23][LR]\s|^[X]\s" \
    -g TrEMOLO/modules/1-FREQUENCY-MULTI-GENERATIONS/test/ref.fasta
```

To see the graph, open `TEST-GRAPH-FRQUENCIES/index.html` file


## Graph

### Generational Frequency Graph

This graph offers a visual representation of TE frequencies across generations. Customize your view by filtering:
* Transposable elements (TEs) of interest,
* Trajectories of frequency evolution — increasing, decreasing, stable, or variable,
* The least number of generations sharing a common position,
* Specific generations you want to focus on,
* The chromosome of interest.

<img src="img/ex1.png">

### Detailed TE Frequency Evolution

By selecting a point on the generational frequency graph, the second graph provides a detailed view of the TE frequency changes across selected generations.

<img src="img/ex2.png">

**Note**: Zero values typically indicate TEs that were not detected in the respective generation.
