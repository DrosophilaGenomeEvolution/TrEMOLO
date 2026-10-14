# TE dans les reads non mappés

Cette branche optionnelle du workflow migré recherche les TE dans les reads
sans aucun alignement sur `DATA.GENOME`. Elle fournit des observations sans
position génomique, indépendantes des appels INSIDER/OUTSIDER, de l'intégration
dans le génome et du calcul de fréquence.

## Activation

```yaml
CHOICE:
  PIPELINE:
    UNMAPPED_TE: true
PARAMS:
  UNMAPPED_TE:
    MIN_PIDENT: 80
    MIN_ALIGNED_BP: 80
    MAX_EVALUE: 1e-10
    FULL_LENGTH_COVERAGE: 80
```

Ces entrées s'ajoutent aux sections existantes du YAML. `DATA.SAMPLE` est
obligatoire. La branche est désactivée si `UNMAPPED_TE` est absent ou `false`.
Elle peut aussi être utilisée avec `OUTSIDER_VARIANT: false` : seul le mapping
nécessaire est alors lancé. La cible explicite `unmapped_te` permet de lancer
ce catalogue seul avec la même commande Snakemake, en remplaçant `all`.

Avec Snakemake 5.10, changer un paramètre ne relance pas automatiquement les
sorties déjà présentes. Après modification des seuils, ajouter
`--forcerun blast_te_against_unmapped_reads annotate_unmapped_te` à la commande
et placer la cible `all` avant cette option. Cela recalcule le BLAST, les tables
et les étapes qui en dépendent. Un changement des options de mapping nécessite
également de relancer le mapping et l'extraction des reads.

## Sélection et alignement

La sélection utilise le SAM initial de minimap2 : le BAM OUTSIDER final retire
déjà les reads non mappés (`-F 4`). Un read est retenu si un enregistrement porte
le flag unmapped et qu'aucun enregistrement du même identifiant n'est mappé,
y compris les alignements secondaires et supplémentaires. Les reads partiellement
mappés et les alignements de faible qualité restent exclus de ce catalogue.
Les options de mapping restent celles de `PARAMS.OUTSIDER_VARIANT.MINIMAP2`.
Ne pas utiliser une option minimap2 qui supprime les enregistrements non mappés
(par exemple `--sam-hit-only`) : ils ne seraient plus disponibles pour la sélection.

Les séquences complètes sont reprises du FASTQ original, dans son orientation,
en conservant les identifiants dans une table de correspondance. Les noms du FASTA
sont des identifiants techniques pour éviter les ambiguïtés d'en-tête BLAST.
Les identifiants sélectionnés doivent être uniques dans le FASTQ.

L'extraction utilise un index d'identifiants en mémoire avec un budget estimé de
256 MiB. Si ce budget est dépassé, elle bascule sur une base SQLite temporaire
avec opérations par lots, sans conserver tous les noms en mémoire. Les lots SQL
contiennent au plus 900 noms ; les lots FASTQ du repli sur disque sont limités à
1 000 reads ou environ
4 MiB de données Python, hors cas d'un read individuel plus grand. Dans le cas
en mémoire, les reads sont traités directement, sans copier les séquences des
reads non candidats. Le CLI `extract --name-memory-mb N` permet d'ajuster ce budget ; `0` force
le traitement sur disque.

Deux parcours du SAM restent nécessaires pour exclure tout read ayant aussi un
alignement mappé, quel que soit l'ordre des enregistrements, sans stocker tous les
noms mappés. Le second parcours est omis si aucun candidat non mappé n'existe.
Le FASTQ est parcouru une fois ; les séquences, leur ordre, les identifiants de
sortie et les statistiques restent identiques à l'extraction précédente.
Cette étape reste mono-thread ; l'optimisation réduit les recherches SQL et les
copies de séquences. Le BLAST suivant utilise les threads du workflow.

La banque TE est alignée par `blastn -task blastn` sur une base des reads non
mappés. Le nombre maximal de sujets est fixé au nombre de reads sélectionnés,
pour éviter la limite implicite de 500 sujets. Les HSP sont regroupés avec le
même algorithme que l'annotation des TE résidents : famille, orientation,
ordre et proximité des fragments. Des copies répétées peuvent produire
plusieurs matches dans un read. Les familles alternatives sont conservées.
La couverture du consensus utilise l'union des intervalles, sans compter deux
fois les portions chevauchantes. L'identité est une moyenne pondérée par la
longueur des HSP retenus ; les HSP redondants sont signalés séparément.

## Fichiers dans `UNMAPPED_TE/`

- `unmapped_reads.fasta` : tous les reads sélectionnés, même sans TE détecté.
- `unmapped_reads.tsv` : `query_id` technique, `read_id` original, `read_length`.
- `UNMAPPED_TE.tsv` : un match TE retenu par ligne, alternatives comprises.
- `UNMAPPED_TE_HSP.tsv` : tous les HSP retournés par BLAST, avec les rejets locaux
  et les redondances ; les hits non retournés par BLAST ne sont pas représentés.
- `extraction.json` : effectifs de sélection et exclusions pour alignement mappé.
- `summary.json` : effectifs, ambiguïtés, paramètres et convention de coordonnées.
- `README.md` : copie de cette description.
- `TE_vs_reads.bln`, `BLAST_DB/` : résultats BLAST bruts et base des reads.

Les deux tables TE ont une ligne d'en-tête et les noms originaux de la banque.
Toutes les coordonnées sont **0-based, fin exclue**, dans le read original
ou dans le consensus TE. Aucun champ ne représente une position génomique.
En l'absence de reads non mappés, BLAST est ignoré et les tables ont uniquement
leurs en-têtes ; en l'absence de match retenu, le catalogue reste également vide.

## Colonnes de `UNMAPPED_TE.tsv`

| Colonne | Description |
|---|---|
| `match_id` | Identifiant du match regroupant les fragments compatibles. |
| `component_id` | Groupe de matches fortement chevauchants dans le read. Ce n'est pas un locus génomique. |
| `read_id` | Identifiant original, premier champ de l'en-tête FASTQ. |
| `read_length` | Longueur complète du read en bases. |
| `te_family` | Nom original du consensus de la banque TE. |
| `read_start`, `read_end` | Bornes du match dans le read, incluant les espaces entre fragments. |
| `strand` | Orientation du TE par rapport au read original (`+` ou `-`). |
| `te_length` | Longueur du consensus TE. |
| `te_covered_bp` | Bases du consensus couvertes, par union des fragments. |
| `te_coverage` | Couverture du consensus en pourcentage (0–100). |
| `read_aligned_bp` | Bases du read couvertes par union des intervalles des fragments. |
| `identity` | Identité moyenne pondérée en pourcentage. |
| `bitscore` | Somme des bit scores des fragments retenus. |
| `best_evalue` | Plus petite e-value des fragments retenus. |
| `fragment_count` | Nombre de fragments non redondants regroupés. |
| `coverage_status` | `near_full_length` si la couverture atteint le seuil configuré, sinon `partial`. Ne garantit pas un TE fonctionnel ou intact. |
| `assignment` | `primary` pour le représentant du composant, `alternative` pour les autres matches ; ils sont tous conservés. |
| `candidate_count` | Nombre de matches dans le composant, éventuellement de même famille. |
| `read_segments` | Union des intervalles effectivement couverts, sous la forme `début-fin;début-fin`. |

## Colonnes de `UNMAPPED_TE_HSP.tsv`

`fragment_id` identifie un HSP ; `match_id` désigne son match ou `.` s'il n'est
pas associé à un match retenu. `read_id`, `read_length`, `te_family`,
`read_start`, `read_end`, `strand`, `te_length` ont le même sens que ci-dessus.
`te_start`, `te_end` sont les bornes sur le consensus. `aligned_bp` est la
longueur de l'alignement BLAST, `identity` son identité, `evalue` son e-value,
`bitscore` son bit score. `filter_status` vaut `accepted` ou `rejected` ;
`filter_reason` indique `identity`, `aligned_bp`, `evalue`, ou `.` sans rejet.
`redundant_with` donne l'identifiant du fragment représentatif, ou `.`.

## Paramètres

Les paramètres sont indépendants des seuils d'appels INSIDER/OUTSIDER :

| Paramètre | Défaut | Rôle |
|---|---:|---|
| `MIN_PIDENT` | 80 | Identité minimale de chaque HSP, en %. |
| `MIN_ALIGNED_BP` | 80 | Longueur minimale de chaque HSP en bases. |
| `MAX_EVALUE` | 1e-10 | E-value maximale pour BLAST et la sélection des HSP. |
| `FULL_LENGTH_COVERAGE` | 80 | Couverture donnant le statut `near_full_length`, en %. |
| `MIN_CHAIN_GAP` | 200 | Tolérance minimale d'espace entre fragments compatibles, en bases. |
| `MAX_GAP_FRACTION` | 0.2 | Tolérance d'espace proportionnelle à la longueur du consensus. La tolérance appliquée est le maximum des deux valeurs. |
| `MAX_OVERLAP` | 20 | Chevauchement maximal entre fragments chaînés, en bases sur le read et le consensus. Une tolérance limitée évite de réunir des répétitions de courts fragments TE. |
| `AMBIGUITY_OVERLAP` | 0.8 | Chevauchement réciproque minimal pour grouper les matches en composant. |

## Interprétation

Un hit est une observation de similarité avec la banque dans un read non mappé.
Il n'établit pas une nouvelle insertion, sa localisation, une copie indépendante
par read, ni une fréquence allélique. Divergence, séquence absente de l'assemblage,
qualité des reads et contamination peuvent également expliquer ces observations.
Le catalogue annexe n'alimente pas `TE_INFOS.bed`, les fichiers de fréquence ou
l'intégration des TE dans le génome.
