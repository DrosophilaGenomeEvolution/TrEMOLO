# MODULE ANALYSYS TE BLAST

This module enables the visualization of BLAST results concerning the newly detected transposable element insertions. It allows for the visual identification of specific structures such as LTR recombinations, transposable elements (TEs) inserted within other TEs, or more complex structures like clusters of TEs. This tool is crucial for genomic researchers aiming to deeply analyze the dynamics of TE insertions.

## Run Build Data

Execute the module with the following command :

```
singularity exec TrEMOLO.simg TrEMOLO/modules/2-MODULE_TE_BLAST/server/scripts/buildData.sh <work-directory-path>
```

For running tests

```bash
#after running the following command 
# singularity exec TrEMOLO.simg snakemake --snakefile TrEMOLO/run.snk --configfile TrEMOLO/test/tmp_config.yml

# run test
singularity exec TrEMOLO.simg TrEMOLO/modules/2-MODULE_TE_BLAST/server/scripts/buildData.sh work_test
```

## RUN SERVER

```bash
#get dependancies
singularity exec TrEMOLO.simg npm install blessed

#run
singularity exec TrEMOLO.simg bash TrEMOLO/module start 2
```

You can change the PORT numbers in the `config.yaml` file
