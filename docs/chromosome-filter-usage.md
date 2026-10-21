# Emplacements des filtres de chromosomes

Inventaire du checkout local au 21 octobre 2026, après ajout des messages de progression.
Il porte sur `CHROM_KEEP` et les options regex de chromosomes, pas sur toutes les regex de parsing du projet.
Les lignes de configuration, validation, transmission, application et diagnostic sont incluses.
Le workflow actuel est `workflow/Snakefile`. Le workflow historique et les utilitaires indépendants sont identifiés séparément.
Les numéros ci-dessous correspondent à cette révision de travail ; ils peuvent changer après une modification.

**130 lignes dans 18 fichiers.**

## Workflow actuel et moteur commun

### [lib/python/parsing/extract_region_reads_vcf.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/parsing/extract_region_reads_vcf.py:103)

```text
103: from chromosome_filter import compile_chromosome_filter, validate_vcf_chromosomes
105:     chromosome_filter = compile_chromosome_filter(args.chrom)
106:     validate_vcf_chromosomes(args.vcf_file, chromosome_filter)
120: def regexInList(value, liste):
121:     return chromosome_filter.search(value) is not None
131: #regexInList(chrom, chromList)
133: chromList  = args.chrom
136: print("[" + str(sys.argv[0]) + "]", "chromosome liste : ", chromList)
181:             if typeV[0] == "<" and typeV not in typeList and regexInList(chrom, chromList) and expEnd != None and expRName != None :
```

### [lib/python/parsing/get_seq_vcf.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/parsing/get_seq_vcf.py:71)

```text
71: from chromosome_filter import compile_chromosome_filter, validate_vcf_chromosomes
108:     chromosome_filter = compile_chromosome_filter(args.chrom)
109:     validate_vcf_chromosomes(args.vcf_file, chromosome_filter)
124: chrom_list  = args.chrom
131: print("[" + str(sys.argv[0]) + "] : regex contig list : ", chrom_list)
135: def regex_in_list(value, liste):
136:     return chromosome_filter.search(value) is not None
160:             if regex_in_list(contig, chrom_list):
169:         print("[" + str(sys.argv[0]) + "] : ERROR please change your regex (-c option) > regex=" + str(args.chrom) + " :: list=" + str(chrom_list))
199:                     if regex_in_list(chrom, chrom_list) and condition :
240:                             if regex_in_list(chrom, chrom_list) and condition :
279:                     if regex_in_list(chrom, chrom_list) and condition :
```

### [lib/python/parsing/sniffles_vcf.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/parsing/sniffles_vcf.py:11)

```text
11: from chromosome_filter import compile_chromosome_filter, validate_vcf_chromosomes
72: def selected(record, args):
75:             compile_chromosome_filter(args.chrom).search(record["chrom"]) is not None)
79:     validate_vcf_chromosomes(args.vcf_file, compile_chromosome_filter(args.chrom))
84:             if not selected(record, args):
104:     validate_vcf_chromosomes(args.vcf_file, compile_chromosome_filter(args.chrom))
112:         if not selected(record, args) or (ids is not None and record["id"] not in ids):
```

### [lib/python/workflow/annotate_resident_te.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/workflow/annotate_resident_te.py:18)

```text
18: from chromosome_filter import compile_chromosome_filter
251:     chrom_pattern: re.Pattern,
297:             if chrom_pattern.search(values["sseqid"]) is None:
891:     parser.add_argument("--chrom-regex", default=".")
924:         chrom_pattern = compile_chromosome_filter(args.chrom_regex)
926:         parser.error(f"invalid chrom-regex: {error}")
932:         chrom_pattern,
```

### [lib/python/workflow/chromosome_filter.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/workflow/chromosome_filter.py:12)

```text
12:         raise ValueError("CHROM_KEEP must be a nonempty regular expression")
40:         raise ValueError("CHROM_KEEP contains an empty comma-separated pattern")
46:         self.patterns = tuple(re.compile(pattern) for pattern in split_patterns(expression))
50:             match = pattern.search(chromosome)
56: def compile_chromosome_filter(expression):
58:         raise ValueError("CHROM_KEEP must be a nonempty regular expression")
67:         raise ValueError(f"Invalid CHROM_KEEP regular expression {expression!r}: {error}") from error
79:     if names and not any(pattern.search(name) for name in names):
82:             "CHROM_KEEP matches no chromosome name in the VCF. "
91:         if fields and pattern.search(fields[0]):
101:         pattern = compile_chromosome_filter(args.regex)
```

### [workflow/Snakefile](/home/mourdas/Documents/5-GIT/TrEMOLO/workflow/Snakefile:16)

```text
16: from chromosome_filter import compile_chromosome_filter
```

### [workflow/rules/outsider_te.smk](/home/mourdas/Documents/5-GIT/TrEMOLO/workflow/rules/outsider_te.smk:89)

```text
89: OUTSIDER_CHROM_KEEP = OUTSIDER_TE_PARAMS.get("CHROM_KEEP", ".")
91:     compile_chromosome_filter(OUTSIDER_CHROM_KEEP)
93:     raise ValueError(f"PARAMS.OUTSIDER_VARIANT.TE_DETECTION.CHROM_KEEP: {error}")
94: CHROMOSOME_FILTER_SCRIPT = str(PIPELINE_ROOT / "lib/python/workflow/chromosome_filter.py")
282:         chromosome_filter=CHROMOSOME_FILTER_SCRIPT,
304:         chrom=OUTSIDER_CHROM_KEEP,
712:         chromosome_filter=CHROMOSOME_FILTER_SCRIPT,
735:         chrom=OUTSIDER_CHROM_KEEP,
762:             bedtools sort -i {input.sniffles:q} | python3 {input.chromosome_filter:q} --regex {params.chrom:q} \
768:             bedtools sort -i {input.direct:q} | python3 {input.chromosome_filter:q} --regex {params.chrom:q} \
784:             bedtools sort -i {input.soft:q} | python3 {input.chromosome_filter:q} --regex {params.chrom:q} \
805:             bedtools sort -i {input.hard:q} | python3 {input.chromosome_filter:q} --regex {params.chrom:q} \
827:                 | python3 {input.chromosome_filter:q} --regex {params.chrom:q} > {output.hard_only:q} || true
835:                 python3 {input.chromosome_filter:q} --regex {params.chrom:q} --input "$evidence" >> {output.merged:q}
```

### [workflow/rules/te_genome.smk](/home/mourdas/Documents/5-GIT/TrEMOLO/workflow/rules/te_genome.smk:43)

```text
43: TE_GENOME_CHROM_REGEX = TE_GENOME_PARAMS.get("CHROM_KEEP", ".")
45:     compile_chromosome_filter(TE_GENOME_CHROM_REGEX)
47:     raise ValueError(f"PARAMS.TE_GENOME.CHROM_KEEP is not a valid regular expression: {error}")
158:         chromosome_filter=CHROMOSOME_FILTER_SCRIPT,
169:         chrom_regex=TE_GENOME_CHROM_REGEX,
193:             --target {wildcards.target:q} --chrom-regex {params.chrom_regex:q} \
```

## Export GFA et autres utilitaires indépendants

### [lib/bash/correction_genome_common.sh](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/bash/correction_genome_common.sh:59)

```text
59:     local chromosome_regex=$7
66:         -v chromosome_regex="$chromosome_regex" '
106:                 if ($1 !~ chromosome_regex) {
108:                     reason="chromosome_filter"
```

### [lib/bash/correction_genome_transloc_diff_chrom.sh](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/bash/correction_genome_transloc_diff_chrom.sh:14)

```text
14: CHROMOSOME_REGEX="."
36:   -r, --chrom_regex REGEX       Only use query chromosomes matching REGEX [.]
66:         -r|--chrom_regex)
68:             CHROMOSOME_REGEX=$2
210:     "$MINIMUM_SIZE" "$MINIMUM_MAPQ" "$MINIMUM_IDENTITY" "$CHROMOSOME_REGEX"
```

### [lib/bash/correction_genome_transloc_sim_chrom.sh](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/bash/correction_genome_transloc_sim_chrom.sh:13)

```text
13: CHROMOSOME_REGEX="."
34:   -r, --chrom_regex REGEX       Only use query chromosomes matching REGEX [.]
60:         -r|--chrom_regex)
62:             CHROMOSOME_REGEX=$2
196:     "$MINIMUM_SIZE" "$MINIMUM_MAPQ" "$MINIMUM_IDENTITY" "$CHROMOSOME_REGEX"
```

### [lib/python/locus/build_te_gfa.py](/home/mourdas/Documents/5-GIT/TrEMOLO/lib/python/locus/build_te_gfa.py:177)

```text
177:           min_identity=80, min_coverage=80, min_aligned_bp=100, selected_tes=None, chrom_regex=None,
180:         chrom_pattern = re.compile(chrom_regex) if chrom_regex is not None else None
196:     if chrom_pattern is not None:
197:         genome = {chrom: sequence for chrom, sequence in genome.items() if chrom_pattern.search(chrom)}
212:                                                   set(genome) if chrom_pattern is not None or only_called_chromosomes else None)
338:         chrom_regex=chrom_regex, chromosomes=list(genome),
362:     parser.add_argument("--chrom-regex", help="Additional Python regex searched on reference sequence names")
365:                                  help="Include reference chromosomes absent from TE_INFOS.bed; still applies --chrom-regex")
385:                      args.min_aligned_bp, selected_tes, args.chrom_regex, args.only_called_chromosomes)
```

## Workflow historique

### [rules.snk](/home/mourdas/Documents/5-GIT/TrEMOLO/rules.snk:27)

```text
27:         CHROM_KEEP       = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
256:         CHROM_KEEP       = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
669:             grep -E `echo "{params.CHROM_KEEP}" | tr "," "|"` {input.genome}.fai | grep -w -f {params.work_directory}/tmp_chrom_with_TE.txt 2> /dev/null | awk 'BEGIN{{print "var chroms = {{"}} {{print "   \\""$1"\\" : "$2","}} END{{print "}}"}}'  > {params.work_directory}/js/chroms.js || fail_msg "ERROR : NO CHROM IN REPORT\\n"
670:             # grep -E `echo "{params.CHROM_KEEP}" | tr "," "|"` {input.genome}.fai | grep -w -f {params.work_directory}/tmp_chrom_with_TE.txt 2> /dev/null | awk 'BEGIN{{print "const all_data = {{"}} {{print "   \\""$1"\\" : data"$1","}} END{{print "}}"}}'   > {params.work_directory}/js/data_ALL_FREQ_CHROM.js || fail_msg "ERROR : NO TE FOR FREQUENCY IN REPORT\\n"
2521:         CHROM_KEEP             = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
3005:         CHROM_KEEP             = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
3036:         chrom_keep=`echo "{params.CHROM_KEEP}" | tr "," "|"`
3038:         echo "CHROM KEEP : $chrom_keep;" | tee {log}.out;
3065:             bedtools sort -i {input.sniffles_TE_bed} | grep -E "$chrom_keep" | bedtools cluster -d 100 | awk 'BEGIN{{TE=""; num=""}} {{split($4, sp, "|"); if(num==$7 && sp[1]==TE){{TE=sp[1]; num=$7;}}else{{print $0; TE=sp[1]; num=$7;}} }}' | cut -f 1-6 \
3072:             bedtools sort -i {input.TrEMOLO_TE_bed} | grep -E "$chrom_keep" | bedtools cluster -d 100 | awk 'BEGIN{{TE=""; num=""}} {{split($4, sp, "|"); if(num==$7 && sp[1]==TE){{TE=sp[1]; num=$7;}}else{{print $0; TE=sp[1]; num=$7;}} }}' | cut -f 1-6 \
3092:         NB_TE_sn=`test -s {params.work_directory}/OUTSIDER/TE_DETECTION/POSITION_START_TE.bed && grep -E "$chrom_keep" {params.work_directory}/OUTSIDER/TE_DETECTION/POSITION_START_TE.bed -c || echo 0`
3093:         NB_TE_Tr=`test -s {params.work_directory}/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.bed && grep -E "$chrom_keep" {params.work_directory}/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.bed -c || echo 0`
3099:         NB_TE_Tr=`test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TrEMOLO_TE_OUTSIDER_CLUSTER_100.bed && grep -E "$chrom_keep" {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TrEMOLO_TE_OUTSIDER_CLUSTER_100.bed -c || echo 0`
3100:         NB_TE_sn=`test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed && grep -E "$chrom_keep" {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed -c || echo 0`
3117:             bedtools sort -i {input.soft_te} | grep -E "$chrom_keep" | bedtools cluster -d 100 | awk 'BEGIN{{TE=""; num=""}} {{split($4, sp, "|"); if(num==$7 && sp[1]==TE){{TE=sp[1]; num=$7;}}else{{print $0; TE=sp[1]; num=$7;}} }}' | cut -f 1-6 > {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_SOFT_CLUSTER_100.bed
3119:         echo "      NB TE SOFT : `grep -E "$chrom_keep" {input.soft_te} | wc -l`" | tee -a {log}.out;
3147:             bedtools sort -i {params.work_directory}/OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.bed | grep -E "$chrom_keep" | bedtools cluster -d 100 | awk 'BEGIN{{TE=""; num=""}} {{split($4, sp, "|"); if(num==$7 && sp[1]==TE){{TE=sp[1]; num=$7;}}else{{print $0; TE=sp[1]; num=$7;}} }}' | cut -f 1-6 > {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_HARD_CLUSTER_100.bed
3150:         echo "      NB TE HARD ; `grep -E "$chrom_keep" {input.hard_te} -c`" | tee -a {log}.out
3151:         echo "      NB TE HARD CLUSTERED (-w 100) ; `grep -E "$chrom_keep" {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_HARD_CLUSTER_100.bed -c`" | tee -a {log}.out
3179:             grep -w -v -f {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_ID_HARD_commun.txt {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_HARD_CLUSTER_100.bed | grep -E "$chrom_keep" > {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_HARD_NOT_FOUND_IN_SOFT_sniffles_and_TrEMOLO.bed
3184:         test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed && cat {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed | grep -E "$chrom_keep" >> {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed
3185:         test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_TrEMOLO_NOT_FOUND_IN_sniffles.bed && cat {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_TrEMOLO_NOT_FOUND_IN_sniffles.bed | grep -E "$chrom_keep" >> {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed
3193:         test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_SOFT_NOT_FOUND_IN_sniffles_and_TrEMOLO.bed && cat {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_SOFT_NOT_FOUND_IN_sniffles_and_TrEMOLO.bed | grep -E "$chrom_keep" >> {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed
3194:         test -s {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_HARD_NOT_FOUND_IN_SOFT_sniffles_and_TrEMOLO.bed && cat {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_HARD_NOT_FOUND_IN_SOFT_sniffles_and_TrEMOLO.bed | grep -E "$chrom_keep" >> {params.work_directory}/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed
3241:         CHROM_KEEP             = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
3378:         CHROM_KEEP             = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
3555:         CHROM_KEEP             = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
3591:         python3 ${{path_to_pipline}}/lib/python/parsing/get_seq_vcf.py {params.GET_SEQ_REPORT_OPTION} -c '{params.CHROM_KEEP}' {input.vcf} {output.snif_seqs} 2>> {log.DETECTION}.err ;
4230:         CHROM_KEEP     = config["PARAMS"]["OUTSIDER_VARIANT"]["TE_DETECTION"]["CHROM_KEEP"],
4267:         chrom_keep=`echo "{params.CHROM_KEEP}" | tr "," "|"`
4271:         cut -f 1,7,8 {params.work_directory}/INSIDER/TE_DETECTION/ALL_TE_1.bln | grep -E "$chrom_keep" > {params.work_directory}/INSIDER/TE_DETECTION/ALL_TE_1.bed
```

## Configurations présentes dans ce dépôt

### [config.yaml](/home/mourdas/Documents/5-GIT/TrEMOLO/config.yaml:54)

```text
54:         CHROM_KEEP: "."          # Python regex on chromosome name only; exact dm6: '^([23][RL]|[X4Y])$'
83:             CHROM_KEEP: "." # Python regex on chromosome name only; exact dm6: '^([23][RL]|[X4Y])$'; '.' keeps all
```

### [config_INSIDER.yaml](/home/mourdas/Documents/5-GIT/TrEMOLO/config_INSIDER.yaml:53)

```text
53:         CHROM_KEEP: "."
79:             CHROM_KEEP: "." # regular expresion of chromosome; exemple  for Drosophila  "[23][RL],4,X" ; Put "." for keep all chromosome
```

### [config_OUTSIDER.yaml](/home/mourdas/Documents/5-GIT/TrEMOLO/config_OUTSIDER.yaml:53)

```text
53:         CHROM_KEEP: "."
82:             CHROM_KEEP: "." # regular expresion of chromosome; exemple  for Drosophila  "[23][RL],4,X" ; Put "." for keep all chromosome
```

### [config_dm6_G73_migrated.yaml](/home/mourdas/Documents/5-GIT/TrEMOLO/config_dm6_G73_migrated.yaml:64)

```text
64:       CHROM_KEEP: "^[23][RL],^[X4Y]"
```

### [tests/workflow/refactor_config.yml](/home/mourdas/Documents/5-GIT/TrEMOLO/tests/workflow/refactor_config.yml:51)

```text
51:     CHROM_KEEP: "."
78:       CHROM_KEEP: "^[23][RL]_,^[X4Y]_"
```
