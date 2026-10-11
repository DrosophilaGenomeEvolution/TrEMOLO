# 0023 — Catalogue annexe des TE dans les reads non mappés

## Statut

Accepté.

## Contexte

Les reads sans alignement génomique peuvent contenir des séquences TE qui ne
figurent pas dans les appels localisés. Le BAM OUTSIDER les exclut déjà.
Une observation dans ces reads ne suffit pas à localiser une insertion ou à
calculer sa fréquence dans la population.

## Décision

Ajouter `CHOICE.PIPELINE.UNMAPPED_TE`, désactivé par défaut. Utiliser le SAM
initial pour sélectionner les identifiants unmapped, puis exclure tout identifiant
qui possède un autre alignement mappé. Récupérer la séquence complète et son
orientation depuis le FASTQ original. Les identifiants de sélection sont stockés
temporairement dans SQLite pour éviter un ensemble de tous les reads en mémoire.

Aligner la banque TE sur les reads sélectionnés par BLAST, sans limite implicite
de 500 sujets. Réutiliser les fonctions de normalisation, de regroupement des
fragments et d'ambiguïté de l'annotation des TE résidents, avec des paramètres
propres aux reads. Conserver chaque famille alternative dans la table finale,
les HSP rejetés dans une table de preuve et les reads sans TE dans le FASTA.

Publier des tables avec en-têtes, noms TE originaux et coordonnées sur le read
0-based, fin exclue. Les fichiers ne sont pas des entrées de `TE_INFOS.bed`, du
calcul de fréquence ou de l'intégration génomique. Les durées des trois nouvelles
étapes participent au graphique existant du REPORT si cette branche est activée.
La production de GFA est différée et ne fait pas partie de ce changement.

## Validation

Les tests doivent couvrir les alignements secondaires/supplémentaires, les
séquences originales compressées, les fragments, les répétitions, les familles
ambiguës, les orientations, les noms TE originaux et les catalogues vides.
Un test du workflow autonome doit retrouver des TE connus dans des reads
synthétiques non mappés. Le workflow complet doit conserver les résultats
localisés, les fréquences et les génomes intégrés de la validation précédente.
