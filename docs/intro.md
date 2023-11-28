# Hydra-genetics reports module

The reports module consists of a set of rules for creating various reports of results within the context of Hydra Genetics. Reports included at the moment are:

- [CNV report](/reports/#cnvs)

## How to use the module

Include the module in your Hydra Genetics workflow by loading the module in your Snakefile:

```python
from hydra_genetics.utils.misc import get_module_snakefile

module reports:
    snakefile: get_module_snakefile(config, "hydra-genetics/reports", path="workflow/Snakefile", tag="v0.2.0")
    config: config

use rule * from reports as reports_*
```
