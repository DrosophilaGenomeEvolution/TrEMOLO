# Minimal analysis configurations

Copy one YAML, replace all `/absolute/path/...` values and select a new output
path. These examples use the current scientific image and leave optional
resident annotation and integration disabled. OUTSIDER and combined examples
enable `UNMAPPED_TE: true` to search unmapped reads for TE matches. The
assembly-only INSIDER example keeps it false because it has no `DATA.SAMPLE`;
add reads before enabling that option.

| File | Required inputs |
| --- | --- |
| [outsider.yaml](outsider.yaml) | Mapping assembly, reads, original TE database |
| [insider.yaml](insider.yaml) | Assembly, reference assembly, original TE database |
| [combined.yaml](combined.yaml) | Assembly, reference assembly, reads, original TE database |

From the parent of the checkout:

```bash
cp TrEMOLO/examples/config/outsider.yaml my-analysis.yaml
# Edit my-analysis.yaml, then preview the jobs:
singularity exec TrEMOLO-3.simg snakemake \
  --snakefile TrEMOLO/workflow/Snakefile \
  --configfile my-analysis.yaml --cores 8 --dry-run all
# Remove --dry-run to execute.
```

Reads should correspond to the input assembly/sample. Use compatible chromosome
names and reference coordinates when combining analyses. YAML paths resolve
from the command working directory. For more settings, see the
[configuration reference](../../docs/configuration.md).
