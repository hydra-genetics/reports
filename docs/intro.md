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

### Overriding rules

Sometimes it might be required to overload only a certain input for a rule. This is however currently not possible with Snakemake&mdash;every input has to be redefined when overloading the rule. You can however make use of the `workflow` object.

Given the example above, let's say we want to overload the `tc_file` input and the associated `tc` and `tc_method` parameters for the rule `cnv_html_report`, but use the default for everything else, we could do this:

```python
use rule cnv_html_report from reports as reports_cnv_html_report with:
    input:
        json=workflow.get_rule("reports_cnv_html_report").input.json,
        html_template=workflow.get_rule("reports_cnv_html_report").input.html_template,
        css_files=workflow.get_rule("reports_cnv_html_report").input.css_files,
        js_files=workflow.get_rule("reports_cnv_html_report").input.js_files,
        tc_file="{sample}_{type}.tc.txt",
    params:
        include_table=workflow.get_rule("reports_cnv_html_report").params.include_table,
        include_cytobands=workflow.get_rule("reports_cnv_html_report").params.include_cytobands,
        tc=lambda wc: open("{sample}_{type}.tc.txt".format(**wc)).read().strip(),
        tc_method="override",
```

If wildcards are being used, make sure that those wildcards are also defined in the output section of that particular rule for things to function as expected.

!!! note
    If the default rule arguments point to files within the module itself, then Snakemake will try to access those from its source cache. If you are running Snakemake with Apptainer/Singularity the path to the cache must be bound for this to work. [The location is determined by the appdirs package](https://snakemake.readthedocs.io/en/stable/executing/cli.html#important-environment-variables), and on a unix/linux system this is determined by the environment variable `$XDG_CACHE_HOME`. The exact path can be obtained by finding the `user_cache_dir` entry after running `python -m appdirs`.

## Versioning caveat

The module makes use of the release-please Github action, meaning that the templates get tagged with the version number every time a new release is published. This means that only reports generated from the main branch are guaranteed to be 100% accurate when it comes the version. Any commit between releases will have the version of the most recent release associated with it.
