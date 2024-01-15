# MODULE SCATTER FREQUENCY TE TrEMOLO

## PREPARE INPUT

pour fonctionné, ce module à besoin d'un fichier qui contient les chemins des différentes analyses effectué avec TrEMOLO (minimum +v2.5.1)

exemple : 

```
work_diretory_1
/path/to/work_diretory_2
work_diretory_3/
```

l'ordre des chemins est important elle désigne l'ordre par lequel le module considérera comme ayant des géneration de la plus ancienne à la plus récente.

Une autre solutions vous permet de spécifier vous même le temps des génerations (format work_directory:G[NUMBER]), exemple :

```
work_diretory_1:G1
/path/to/work_diretory_2:G10
work_diretory_3:G3
```

ce format vous permet d'indiqué l'ordre et le temps de décalge entre les génerations, dans l'exemple précedent la G1 et la géneration la plus ancienne tandis que la G10 la plus récente.

***L'ordre des génerations vous permet au module d'indiquée quelles sont les TE qui augmente, diminue, ou varie au cours des génerations.***


## RUN BUILD GRAPH

```
singularity exec TrEMOLO.simg TrEMOLO/modules/1-FREQUENCY-MULTI-GENERATIONS/buildFrequencyGenerations.sh -i <input-init-file> [-o OUTPUT_NAME_DIRECTORY]
```

