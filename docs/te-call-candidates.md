# Familles candidates et provenance des séquences

Le workflow migré écrit deux tables avec en-têtes et noms TE originaux :

- `TE_CALL_CANDIDATES.tsv` : synthèse des familles candidates, seulement pour les
  événements retenus dans `TE_INFOS.bed` qui possèdent plusieurs familles.
- `TE_CALL_CANDIDATE_EVIDENCE.tsv` : HSP BLAST bruts qui soutiennent ces candidats,
  reliés aux séquences analysées et, lorsqu'il est connu, au read d'origine.

Cette mise à jour ne change ni les familles principales, ni les comptes historiques,
ni les fréquences, le TSD ou les génomes intégrés. Elle ne construit aucun consensus
et ne regroupe pas des événements distincts en locus biologique.

## Synthèse : `TE_CALL_CANDIDATES.tsv`

Les 17 colonnes historiques restent dans le même ordre :

| Colonnes | Description |
|---|---|
| `candidate_group_id` | Identifiant stable du groupe de familles d'un événement. |
| `source` | Branche INSIDER ou OUTSIDER. |
| `chrom`, `start`, `end` | Position de l'appel, coordonnées 0-based, fin exclue. |
| `event_id`, `tremolo_id`, `event_type` | Identifiants de l'événement et type de SV. |
| `reported_te` | Famille principale de `TE_INFOS.bed`. |
| `candidate_te` | Famille candidate de cette ligne. |
| `assignment` | `reported_primary` ou `alternative`. |
| `candidate_rank`, `candidate_count` | Rang et nombre de familles candidates. |
| `evidence_count`, `evidence_fraction` | Compte historique et part de preuves attribuée à la famille ; ce ne sont pas des comptes de reads uniques ni des fréquences alléliques. |
| `evidence_channels` | Canaux ayant contribué aux comptes. |
| `ambiguity_type` | Description des relations de séquence, expliquée ci-dessous. |

Six colonnes sont ajoutées. Leurs valeurs portent sur **tout le groupe**, et sont
donc répétées sur ses lignes de familles :

| Colonne | Description |
|---|---|
| `query_count` | Nombre de couples canal / identifiant de séquence analysée. |
| `sequence_count` | Nombre de séquences d'insertion exactement différentes, selon leur SHA256. |
| `read_count` | Nombre d'identifiants de reads connus et distincts. Les origines VCF et inconnues sont exclues. |
| `vcf_query_count` | Nombre de requêtes correspondant à une séquence rapportée dans le VCF. |
| `provenance_status` | `complete` si les séquences et leurs origines sont connues pour tous les canaux/familles ; `partial` si seule une partie est disponible ; `unavailable` sans preuve détaillée. Une origine VCF connue reste complète sans identifiant de read. |
| `relation_types` | Relations observées, séparées par `;`, ou `.`. |

`ambiguity_type` peut valoir :

- `same_sequence_annotation_ambiguity` : des familles différentes ont des HSP
  chevauchants sur la même requête ; une seule séquence peut recevoir plusieurs
  annotations.
- `same_sequence_multiple_regions` : des familles différentes ont des HSP sans
  chevauchement sur la même requête ; cela expose une structure potentiellement
  composite, sans confirmer ses composants.
- `different_sequence_evidence` : les preuves impliquent plusieurs séquences
  d'insertion non identiques, sans les deux relations précédentes.
- `mixed_sequence_evidence` : plusieurs de ces situations coexistent.
- `unresolved_family_candidates` : les preuves ne permettent pas ces distinctions.

Les relations fines sont `same_query_overlapping_families`,
`same_query_disjoint_families`, `different_insertion_sequences`, `different_reads`
et `identical_insertion_multiple_queries`. Le chevauchement est celui des HSP
bruts : même un chevauchement court est signalé, sans décision d'annotation.

Deux reads peuvent contenir une insertion identique : `read_count` vaut alors 2
et `sequence_count` 1. Inversement, une différence de séquence peut provenir d'une
erreur de lecture. Ces valeurs ne prouvent donc pas des allèles distincts ou un SWAP.

## Preuves : `TE_CALL_CANDIDATE_EVIDENCE.tsv`

Une ligne représente un HSP brut, sans nouvelle sélection scientifique. Les
comptes historiques portent sur des associations de séquences/familles et peuvent
différer du nombre de HSP. Les hits faibles restent visibles : leur présence ne
signifie pas qu'ils passent les seuils de classification du TE.

| Colonnes | Description |
|---|---|
| `candidate_group_id`, `event_id` | Raccord à la synthèse. |
| `evidence_channel` | Canal d'origine du BLAST. |
| `query_id` | En-tête exact de la séquence utilisée comme requête BLAST. |
| `read_id` | Read d'origine, ou `.` s'il n'est pas attribuable. |
| `sequence_origin` | `read_sequence`, `vcf_reported` ou `unknown`. |
| `genomic_alignment_flag`, `genomic_alignment_status` | Flag SAM et statut de l'alignement read–génome ayant fourni cette requête, et non du hit BLAST. Ces colonnes sont ajoutées à la fin de la table. |
| `insertion_sha256`, `insertion_length` | Empreinte et longueur de l'insertion seule, sans flancs, ou `.` sans cette information. L'empreinte respecte l'orientation de la séquence extraite ; aucun complément inverse n'est fusionné implicitement. |
| `candidate_te` | Nom original de la famille alignée. |
| `query_start`, `query_end` | Bornes du HSP sur la requête BLAST, 0-based, fin exclue. Cette requête peut contenir des flancs : ces coordonnées ne sont pas celles du FASTQ complet ni celles de l'insertion seule. |
| `te_start`, `te_end` | Bornes sur le consensus TE, 0-based, fin exclue. |
| `strand` | Orientation relative du TE et de la requête. |
| `identity`, `aligned_bp`, `evalue`, `bitscore` | Identité en %, longueur et scores du HSP. |
| `relation_types` | Relations observées dans le groupe. |

## Origine des séquences et rapports

Les extracteurs écrivent `SEQUENCE_PROVENANCE.tsv` près de leurs FASTA, avec
`query_id`, `event_id`, `read_id`, `sequence_origin`, `genomic_alignment_flag`
et `genomic_alignment_status`. Les anciennes tables à quatre colonnes restent
acceptées : leur statut est `unknown` pour un read et `not_applicable` pour une
séquence VCF. Un flag manquant n'est jamais déduit du nom de read ni du MAPQ.
Une insertion CIGAR est associée au read effectivement utilisé. La séquence
rapportée par le VCF est marquée `vcf_reported` avec `read_id=.` : Sniffles 2 peut
produire un consensus et une identité de read ne peut pas être déduite du dernier
alignement parcouru. Les anciens fichiers `RD_NUMBER.txt` ne remplacent pas cette
nouvelle table de provenance.

### Alignement génomique à l'origine de la séquence

Le statut est celui de **l'enregistrement SAM/BAM effectivement utilisé** pour
extraire la séquence, même si plusieurs enregistrements portent le même nom de
read :

- `primary` : ni le bit `0x100` ni le bit `0x800` n'est positionné ;
- `secondary` : bit `0x100` ;
- `supplementary` : bit `0x800` (alignement supplémentaire, distinct d'un secondaire) ;
- `secondary_and_supplementary` : les deux bits sont positionnés ;
- `unmapped` : bit `0x4`, conservé explicitement s'il est rencontré ;
- `unknown` : aucun flag disponible, notamment pour les résultats antérieurs ;
- `not_applicable` : séquence rapportée par le VCF, sans alignement source unique.

Les insertions CIGAR directes et les insertions extraites près des événements
Sniffles gardent le flag de leur enregistrement source. Pour les soft clips,
la provenance suit exactement le représentant gauche/droit déjà choisi, y
compris le premier en cas d'égalité de longueur. Pour les hard clips récupérés
dans le FASTQ, elle garde le flag du candidat génomique qui a demandé cette
extraction. Les choix de représentants, les séquences, les appels TE et les
fréquences ne changent pas. En mode refondu, le flag est ajouté après les huit
colonnes historiques du BED des insertions directes et après les sept colonnes
historiques du BED des fragments HARD, afin de traverser le clustering sans
rechercher un alignement à partir du seul nom de read.

Le REPORT principal ajoute **Source genomic alignments (queries)** : un décompte
des requêtes extraites par statut pour chaque événement OUTSIDER, tous canaux
confondus. Plusieurs statuts peuvent coexister ; le même alignement/read peut
être compté via plusieurs requêtes ou canaux. Ce n'est ni un nombre de reads
uniques, ni un nombre d'alignements indépendants, ni le support de chaque famille
TE affichée. Des requêtes du même événement sans hit TE retenu peuvent aussi
figurer dans ce décompte. Une origine INSIDER assemblage est affichée comme non
applicable. Le module 2 et les preuves détaillées des appels ambigus affichent
le statut et le flag par requête ; `insertion-sequence-provenance.tsv` les conserve.

Sur une nouvelle analyse, le coût ajouté se limite à copier le flag pendant les
lectures/extractions existantes et à écrire/lire les tables de provenance.
Il n'y a ni étape d'alignement/BLAST supplémentaire, ni nouveau parcours du BAM.
Les rapports ne relisent jamais les BAM pour compléter un flag manquant. Aucun
gain ou pourcentage de surcoût sur données réelles n'a été mesuré pour cet ajout.
Pour d'anciens résultats, afficher `unknown` est possible sans recalcul ; obtenir
des flags absents exige de relancer les extracteurs concernés et leurs sorties
dépendantes, ce qui n'est pas une opération gratuite.

Le REPORT principal expose la relation, les comptes de séquences/reads, l'état de
provenance, les HSP consultables et le téléchargement des deux tables. Le module
`2-MODULE_TE_BLAST` expose aussi l'origine de chaque insertion et écrit
`insertion-sequence-provenance.tsv` quand les données sont disponibles.
Les anciennes tables de candidats à 17 colonnes restent lisibles par le rapport,
avec provenance indisponible. Les appels INSIDER n'acquièrent pas de provenance
inventée : leur classifieur conserve encore une seule meilleure famille.

Pour actualiser une analyse déjà calculée, relancer `all` en plaçant cette cible
avant `--forcerun cluster_outsider_alignment_insertions extract_outsider_sniffles_sequences
build_te_call_candidates`. Les étapes dépendantes seront recalculées ; les
extractions ajoutent les tables de provenance sans modifier les séquences.
