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

## Versioning caveat

The module makes use of the release-please Github action, meaning that the templates get tagged with the version number every time a new release is published. This means that only reports generated from the main branch are guaranteed to be 100% accurate when it comes the version. Any commit between releases will have the version of the most recent release associated with it.
