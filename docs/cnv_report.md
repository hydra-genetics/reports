# CNV report

## Input files

Required files are:

- log<sub>2</sub> ratios from each CNV caller
- Segments from each CNV caller

Optional input files are:

- Filtered and unfiltered VCFs that each contain calls from all included callers for generating a table of CNV calls
- BED files with annotations that should be added to plots
- Germline VCF file for displaying VAF in plots
- Cytoband definition for displaying in the chromosome plot

## Output files

- `reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html`

## Configuration

There are a couple of things that can be customised using the config file.

### Results table

The CNV results table contains CNVs that have been called by the pipeline. In order for the table to be included in the final report, `show_table` under [`cnv_html_report`](/softwares/#configuration) has to be `true`. If this is the case, then both `filtered_cnv_vcfs` and `unfiltered_cnv_vcfs` have to be defined under [`merge_cnv_json`](/softwares/#configuration_2).

### Additional tables

Additional tables can be included in the final report by making use of `extra_tables` under [`cnv_html_report`](/softwares/#configuration). A table should be represented by a tsv file, and the first row will be used as a header for the table. The value of `extra_tables` in the config should be an array of objects, and the objects should look like this:

```yaml
extra_tables:
    - name: Extra table
      description: A description of the table
      path: extra_table.tsv
```


`name` is the name of the table, and will be used as a section heading. `description` is a description of the table and will be displayed as a single paragraph, and `path` is the path to the tsv file from which the table should be created. If the table file is completely empty, the execution will fail with an error. If the table is empty, but it contains a header, the table will be presented. It will however have a message clarifying that there is no data in the table, and that this is how it is meant to be. Wildcards are allowed in `path`, as long as the same wildcards are present in the output file name. By default these wildcards are `sample`, `type` and `tc_method`.

### Cytobands

Cytobands can be represented in the chromosome plot. For these to be included, `cytobands` under [`cnv_html_report`](/softwares/#configuration) has to be `true`, and `cytobands` under [`merge_cnv_json`](/softwares/#configuration_2) should point to a file with cytoband definitions. The format of this file should follow the UCSC cytoband schema ([hg19](https://www.genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema), [hg38](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema)). Currently, files for both hg19 and hg38 are included in the [config directory of the repo](https://github.com/hydra-genetics/reports/tree/develop/config).

In terms of customisation, the colours of the bands can be redefined in `cytoband_config` under [`merge_cnv_json`](/softwares/#configuration_2). For reference, this is the default configuration:

```yaml
merge_cnv_json:
    cytoband_config:
        colors:
            gneg: "#e3e3e3"
            gpos25: "#555555"
            gpos50: "#393939"
            gpos75: "#8e8e8e"
            gpos100: "#000000"
            acen: "#963232"
            gvar: "#000000"
            stalk: "#7f7f7f"
```

Only full-length hexadecimal colours (without alpha channel), as shown above, are supported.

### Custom annotations

Custom annotations can be added to the chromosome plot by specifying one or more bed-files in `annotations` under [`merge_cnv_json`](/softwares/#configuration_2). Only the four first columns of the file will be taken into account, and the value in the name column will be displayed in the plot.

## Customising the template

The template used can be found in [`workflow/templates/cnv_html_report`](https://github.com/hydra-genetics/reports/tree/develop/workflow/templates/cnv_html_report). This will be used by default. If you for some reason would like to customise the template, the input files will have to be redefined when importing the module. Below is an example where template files are redefined, while the input data remains the default:

```snakemake
use rule cnv_html_report from reports as reports_cnv_html_report with:
    input:
        json="reports/cnv_html_report/{sample}_{type}.{tc_method}.merged.json",
        html_template="path/to/custom/template.html",
        css_files=["path/to/custom/css/style.css"],
        js_files=["path/to/custom/js/script-1.js", "path/to/custom/js/script-2.js"],
        tc_file=reports.get_tc_file,
```

The only template file that is strictly required is the html file `html_template`. Both `css_files` and `js_files` can be left out if you so wish, but the functionality will be severly limited without any javascript.
