# Software used in the reports module

## general_json_report

Convert input files to JSON for the general report.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__general_html_report__general_json_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__general_html_report__general_json_report#

### :wrench: Configuration

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__general_json_report#

## general_html_report

Generate a general HTML report for a sample.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__general_html_report__general_html_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__general_html_report__general_html_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__general_html_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__general_html_report#

## cnv_html_report

Generate an HTML report for CNVs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__cnv_html_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnv_html_report__cnv_html_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnv_html_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnv_html_report#

## cnv_json

Convert CNV results from a particular CNV caller to JSON that is compatible with the final report.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__cnv_json#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnv_html_report__cnv_json#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnv_json#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnv_json#

## merge_cnv_json

Merge JSON files from multiple CNV callers and add annotations and other sample specific data.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__merge_cnv_json#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnv_html_report__merge_cnv_json#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__merge_cnv_json#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__merge_cnv_json#
