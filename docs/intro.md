# Hydra-genetics reports module

The reports module consists of a set of rules for creating various reports of results within the context of Hydra Genetics. Reports included at the moment are:

- [CNV report](#cnv-report)

## CNV report

### Input files

Required files are:

- log<sub>2</sub> ratios from each CNV caller
- Segments from each CNV caller

Optional input files are:

- Filtered and unfiltered VCFs that each contain calls from all included callers for generating a table of CNV calls
- BED files with annotations that should be added to plots
- Germline VCF file for displaying VAF in plots

### Output files

- `reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html`
