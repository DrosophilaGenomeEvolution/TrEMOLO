# Tables de fréquence OUTSIDER

Les fichiers `OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv` et
`OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv` possèdent le même en-tête
TSV de **11 colonnes**. Les copies destinées à la lecture sont dans
`PUBLIC/OUTSIDER/FREQUENCY/`, avec les noms originaux des TE. Les fichiers
sources utilisent les identifiants préparés pour les calculs.

| Colonne | Description |
| --- | --- |
| `position` | Locus `chromosome:start-end` du candidat utilisé pour compter les reads sur le génome analysé. Ce n'est pas la position affinée par la recherche de TSD. |
| `event_id` | Identifiant de l'événement, par exemple `TrEMOLO.INS.84` ou `sniffles.INS.25`. |
| `te_family` | Nom du TE attribué au candidat. |
| `support_reads` | Pour une INS : support retenu sans les preuves de clipping, après les corrections du compteur. Pour une DEL : nombre retenu de reads soutenant la délétion, après application du support minimal du caller. |
| `clipped_reads` | Nombre de reads dont la preuve retenue est un soft clip ou un hard clip. Pour une DEL : `.` (non applicable). |
| `total_reads_without_clipped` | Pour une INS : dénominateur de la fréquence sans clipping, après les corrections du compteur. Pour une DEL : total des reads retenus. |
| `te_present_reads` | Pour une INS : support retenu avec clipping, au moins égal au support du caller. Pour une DEL : reads conservant le TE, soit `total_reads - support_reads`. |
| `total_reads` | Dénominateur de la fréquence avec clipping, après les corrections du compteur. |
| `frequency_without_clipped` | Pour une INS : `100 × support_reads / total_reads_without_clipped`. Pour une DEL : `100 × te_present_reads / total_reads`. |
| `frequency_with_clipped` | `100 × te_present_reads / total_reads`. Pour une DEL, identique à la fréquence précédente. |
| `sv_type` | Type `INS` ou `DEL`. Les événements uniquement `HARD` ou `SOFT` n'ont pas de ligne dans ces tables. |

Les fréquences sont exprimées en **pourcentage**, pas entre 0 et 1. Pour une
DEL, elles décrivent la **présence du TE**, pas la fréquence de sa délétion.

## Différence entre les deux fichiers

`FREQUENCY_TE_INS.tsv` résume les preuves par événement et par read. Un read
n'est compté qu'une fois par événement. Le calcul historique conserve la
première preuve dans l'ordre de tri ; une preuve `E` (absence d'insertion)
peut donc masquer une preuve d'insertion ou de clipping pour le même read.
Les comptes peuvent être ajustés au support minimal du caller et ne sont pas
nécessairement des comptages bruts.

`FREQUENCY_TE_INS_PRECISE.tsv` corrige le support des insertions à partir de
`COUNT_TE_IN_RS.txt`. Il contient uniquement les candidats ayant un support
strictement positif. Il remplace `support_reads` par ce support et calcule
`te_present_reads = support_reads + clipped_reads`. Si le support vaut le
dénominateur sans clipping plus un, les deux dénominateurs augmentent de un.
« PRECISE » désigne cette correction de support, pas une précision des
coordonnées. Le calcul peut produire une valeur supérieure à 100 ;
`TE_INFOS.bed` utilise alors les fréquences du fichier non corrigé.

Désactiver `CLIPPED_READS` désactive la détection des appels HARD/SOFT séparés.
Le compteur de fréquence continue à inspecter les soft/hard clips dans les
CIGAR des reads couvrant les candidats INS/DEL.

## Lecture et compatibilité

La première ligne est toujours l'en-tête, même si aucun événement n'est retenu.
Un fichier avec seulement l'en-tête est donc une table vide. Utiliser un lecteur
TSV qui gère les guillemets, par exemple :

```python
import csv
with open("FREQUENCY_TE_INS.tsv", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
```

Pour les noms restaurés contenant `|`, `:`, des espaces ou des guillemets,
les copies publiques utilisent le quoting TSV standard. Le lecteur du workflow
reconnaît et ignore l'en-tête et accepte aussi les anciennes tables sans en-tête.
Les colonnes restent dans leur ordre historique : l'ajout de l'en-tête ne
change ni les comptes, ni les fréquences, ni les identifiants calculés.

`SV_SIZE.tsv`, situé dans le même répertoire, est un intermédiaire de tailles
à deux colonnes (`qseqid`, taille en bases), sans en-tête ; ce n'est pas une
table de fréquence.


## Vue fréquence × profondeur dans le rapport

Le rapport principal ajoute un nuage de points **Frequency versus local depth** :
la profondeur en abscisse et la fréquence affichée en ordonnée. Les filtres du
rapport s'appliquent au graphique. L'axe de profondeur peut être linéaire ou
logarithmique. Deux repères ajustables (initialement 10 reads et 50 %) comptent
les appels dans les quatre combinaisons fréquence/profondeur élevée ou faible.
Ces repères sont des aides exploratoires, pas des seuils de confiance validés.

Chaque point donne au survol le TE, l'événement, les coordonnées, le support du
TE, le dénominateur et la méthode. La table des appels et son export TSV
conservent aussi ces informations. Le graphique fréquence/position existant
inclut désormais la profondeur dans ses infobulles.

- **OUTSIDER** : le rapport sélectionne la même ligne PRECISE, ou la même ligne
  brute de repli, que `TE_INFOS.bed`. Il utilise `total_reads` pour une fréquence
  avec clipping et `total_reads_without_clipped` pour une fréquence sans clipping.
  Le support est celui de la présence du TE ; pour une DEL, il ne s'agit donc
  pas du nombre de reads soutenant la délétion. Ces comptes peuvent être ajustés
  par les règles historiques et ne sont pas une profondeur BAM brute.
- **INSIDER** : `total_depth` de `DEPTH_TE_INSIDER.csv` est la moyenne entière des
  profondeurs échantillonnées aux deux flancs. Le support est `read_support`.
  La correspondance utilise les noms de TE restaurés et les identifiants des
  événements, sans utiliser les positions affinées par TSD pour la jointure.

Une fréquence sans profondeur correspondante reste disponible dans les autres
vues. Elle n'est pas placée artificiellement à profondeur zéro. Les estimations
à profondeur nulle, les comptes invalides, les fréquences hors de [0, 100] et
les données dont la fréquence ou le chromosome ne correspondent pas sont exclus
du nuage et explicités dans `frequency_depth_status`. Aucun dénominateur n'est
reconstruit à partir de la fréquence. Les fréquences et fichiers scientifiques
ne sont pas modifiés par cette vue.

Cette profondeur décrit le support local de l'estimation. Elle ne mesure pas
un enrichissement relatif à une profondeur génomique de référence. Les méthodes
INSIDER et OUTSIDER diffèrent ; utiliser le filtre Source pour les examiner
séparément. Le champ `Coverage` déjà présent dans la table décrit l'alignement
au consensus TE et ne doit pas être confondu avec la profondeur de reads.

Pour régénérer cette vue sur une analyse existante, lancer la cible `report`
avec le même fichier de configuration et le workflow mis à jour. Les tables de
fréquence sont des dépendances déclarées du rapport ; aucun nouveau comptage BAM
n'est nécessaire si elles existent déjà. En appel direct au générateur Python,
fournir `--outsider-frequency`, `--outsider-frequency-precise`,
`--insider-frequency` et `--te-headers` pour les branches disponibles.
