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

It is recommended that only tagged releases of the module are used, and the reason for this is explained in the section [Template version matching](#template-version-matching) below.

### Copying the template

The template for the reports module is located in the [`config/cnv_report_template`](https://github.com/hydra-genetics/reports/tree/main/config/cnv_report_template) directory. This will not be included when loading the model as described above&mdash;it will have to be copied to the main workflow manually. This can be achieved in many different ways, and one way is to use the helper tool <https://download-directory.github.io> and pass in the Github link to the report template. For example, if you want to download v0.2.0 of the module, the final URL would be <https://download-directory.github.io?url=https://github.com/hydra-genetics/reports/tree/v0.2.0/config/cnv_report_template>

!!! important
    Make sure that the version of the template matches with the version of the module that is being imported. More on that below.

## Template version matching

When using this module as a part of a larger workflow, the templates that are used need to be copied to that workflow since they will not be picked up by the module system. To prevent that the module itself and the template goes out of sync, there is a check in place to make sure that they are in fact based on the same version. If not, the workflow execution will fail with an error message.

A caveat with this is that are cases where there could be a mismatch even though the version numbers are identical, for example, if the module is loaded from a commit that lies between releases while the template belongs to the release immediately prior to that same commit. This is why it is recommended that specific releases are used. If a commit between releases is needed for some reason, make sure to copy the template for that same commit.

This is a feature that was added in v0.3.0, and the check is not done in versions prior to this.
