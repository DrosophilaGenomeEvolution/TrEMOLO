# Fenêtres et distances du workflow migré

Les fenêtres utilisées par les commandes BEDTools et les scripts sont
configurables dans `PARAMS.WINDOWS`. Toutes les valeurs sont des entiers
**positifs ou nuls, en bases**. Une option absente conserve sa valeur historique.
Les noms inconnus, valeurs négatives, chaînes et booléens sont refusés au démarrage.

Exemple à intégrer dans le bloc `PARAMS` existant :

```yaml
PARAMS:
  THREADS: 20
  WINDOWS:
    OUTSIDER_INSERTION_CLUSTER_DISTANCE: 50
    OUTSIDER_SOFT_CLIP_WINDOW: 50
    OUTSIDER_HARD_CLIP_CLUSTER_DISTANCE: 100
    OUTSIDER_EVIDENCE_CLUSTER_DISTANCE: 100
    OUTSIDER_EVIDENCE_MATCH_WINDOW: 100
    OUTSIDER_SNIFFLES_READ_MATCH_WINDOW: 30
    OUTSIDER_SNIFFLES_READ_FETCH_WINDOW: 1
    OUTSIDER_SNIFFLES_READ_FLANK_SIZE: 30
    OUTSIDER_CLIPPED_SEQUENCE_PADDING: 100
    OUTSIDER_FREQUENCY_WINDOW: 50
    OUTSIDER_FREQUENCY_REGION_MERGE_DISTANCE: 1
    INSIDER_FREQUENCY_OUTER_MARGIN: 100
    INSIDER_FREQUENCY_INNER_MARGIN: 5
    INSIDER_FREQUENCY_DEPTH_MARGIN: 30
    INSIDER_FREQUENCY_SIZE_WINDOW: 30
    INSIDER_FREQUENCY_BREAKPOINT_DISTANCE: 30
    TE_INFOS_DELETION_EXCLUSION_WINDOW: 30
```

Il suffit de déclarer les options que l'on veut modifier. Les fenêtres ont des
rôles différents ; elles n'ont pas besoin d'avoir la même valeur.

| Option | Défaut | Traitement concerné |
| --- | ---: | --- |
| `OUTSIDER_INSERTION_CLUSTER_DISTANCE` | 50 | `bedtools cluster -d` des insertions extraites des CIGAR. La distance est mesurée entre intervalles BED. |
| `OUTSIDER_SOFT_CLIP_WINDOW` | 50 | Regroupement des soft clips dans `find_all_type_ins.py`, autour du point d'ancrage du groupe ; distance inclusive. |
| `OUTSIDER_HARD_CLIP_CLUSTER_DISTANCE` | 100 | `bedtools cluster -d` des hard clips avant classification TE. |
| `OUTSIDER_EVIDENCE_CLUSTER_DISTANCE` | 100 | Regroupement BEDTools des événements dans chacune des sources Sniffles, CIGAR, SOFT et HARD lors de la fusion. |
| `OUTSIDER_EVIDENCE_MATCH_WINDOW` | 100 | `bedtools window -w` entre ces sources pour repérer les preuves proches. |
| `OUTSIDER_SNIFFLES_READ_MATCH_WINDOW` | 30 | Distance **strictement inférieure** à cette valeur entre l'insertion CIGAR et le point de cassure Sniffles, dans `parse_bam_found_ins.py`. |
| `OUTSIDER_SNIFFLES_READ_FETCH_WINDOW` | 1 | Intervalle BAM `[point − fenêtre, point + fenêtre)` utilisé pour récupérer les reads Sniffles. Ce sont les alignements chevauchant cet intervalle qui sont examinés. Les bornes sont limitées au chromosome. |
| `OUTSIDER_SNIFFLES_READ_FLANK_SIZE` | 30 | Bases de part et d'autre de l'insertion CIGAR extraites avec la séquence du read Sniffles pour classification et préparation TSD. |
| `OUTSIDER_CLIPPED_SEQUENCE_PADDING` | 100 | Bases ajoutées à la longueur du plus grand TE de la banque pour limiter les séquences SOFT/HARD extraites. |
| `OUTSIDER_FREQUENCY_WINDOW` | 50 | Marge des régions utilisées pour extraire le BAM de fréquence **et** fenêtre de recherche des preuves dans `count_outsider_frequency.py`. La valeur est partagée pour ne pas éliminer en amont les reads utiles au compteur. |
| `OUTSIDER_FREQUENCY_REGION_MERGE_DISTANCE` | 1 | `bedtools merge -d` des régions de comptage OUTSIDER. |
| `INSIDER_FREQUENCY_OUTER_MARGIN` | 100 | Extension des frontières INSIDER pour la première sélection des reads. |
| `INSIDER_FREQUENCY_INNER_MARGIN` | 5 | Marge intérieure utilisée pour sélectionner les reads traversant l'insertion. |
| `INSIDER_FREQUENCY_DEPTH_MARGIN` | 30 | Positions de mesure de profondeur à l'extérieur de l'insertion. La même valeur est utilisée dans la préparation des intervalles et le calcul de fréquence. |
| `INSIDER_FREQUENCY_SIZE_WINDOW` | 30 | Tolérance autour de la taille de l'insertion INSIDER pour reconnaître une délétion CIGAR correspondant au site vide ; bornes de taille strictes. |
| `INSIDER_FREQUENCY_BREAKPOINT_DISTANCE` | 30 | Distance stricte d'une extrémité de délétion CIGAR à une frontière INSIDER pour compter le site vide. |
| `TE_INFOS_DELETION_EXCLUSION_WINDOW` | 30 | Marge utilisée pour exclure de `TE_INFOS.bed` un appel OUTSIDER proche d'une délétion INSIDER de la même famille. Les comparaisons d'intersection restent strictes. |

Une valeur nulle ne désactive pas uniformément un traitement : avec les
comparaisons strictes, elle peut empêcher tout appariement ; avec BEDTools,
elle garde la sémantique d'intersection ou de proximité propre à la commande.
Le regroupement `bedtools cluster` peut chaîner des événements successifs :
la portée totale d'un groupe peut donc dépasser `-d`.

`REPORT.LOCUS_WINDOW` reste la portée maximale des groupes de proximité
**affichés dans le rapport**, indépendante des fenêtres d'analyse. De même,
`TSD.SIZE_FLANK`, `LIFT_OFF.FLANK_SIZE`, `LIFT_OFF.MAX_GAP` et les paramètres
de chaînage `TE_GENOME` conservent leurs emplacements dans le YAML.
`GET_SEQ_REPORT_OPTION: "-m 100"` reste un seuil de longueur des séquences
extraites du VCF, pas une fenêtre de proximité.

Pour comparer plusieurs valeurs, utiliser un répertoire de sortie distinct
pour chaque configuration. Sur Snakemake 5, changer les paramètres ne garantit
pas à lui seul le recalcul de résultats existants ; un nouveau
`DATA.WORK_DIRECTORY` permet de repartir sans résultats antérieurs.
