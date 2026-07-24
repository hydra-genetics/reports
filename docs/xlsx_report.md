# XLSX report

The XLSX report compiles SNV, structural variant (SV) and CNV calls for a sample into a single Excel workbook, so that
results can be reviewed without opening several VCFs individually.

## Input files

Required files are:

- A VEP-annotated, phased SNV/indel VCF (plain text or gzipped)
- An SVDB-merged structural variant VCF
- A filtered, annotated CNVkit CNV VCF

## Output files

- `reports/xlsx_reports/{sample}_{type}_combined_report.xlsx`

## Sheets

The workbook contains the following sheets:

Sheet               | Source                     | Description
--------------------|----------------------------|-------------
`SNV`               | SNV VCF                    | SNVs/indels passing filters, excluding TP53 and, if `genes` is configured, restricted to the configured gene list.
`TP53`              | SNV VCF                    | SNVs/indels in the *TP53* gene. Never affected by the `genes` filter.
`Translocations`    | SV VCF                     | `BND` calls on chr4 or chr14.
`SV`                | SV VCF                     | Non-`BND` calls on chr14 that pass filters and have `\|SVLEN\| >= sv_min_len`.
`CNV`               | CNV VCF                    | All CNV calls, with `POS`, `END` and `SVLEN` converted to megabases.
`Software Versions` | `compile_xlsx_report.software_versions` | Tool name and version, parsed from the container tag of each configured tool.

## Configuration

### Filters

The `filters` yaml file, configured under [`compile_xlsx_report`](/softwares/#compile_xlsx_report), controls how the SNV
VCF is parsed and which columns end up in each sheet. It should define:

Key                 | Description
--------------------|-------------
`format_fields`     | FORMAT/sample keys to extract from the SNV VCF (e.g. `GT`, `AF`, `DP`). Missing fields are left blank rather than causing an error.
`vep_info_fields`   | VEP `CSQ` subfields to extract, in the order they appear in the VCF header's `CSQ` definition. Only the first annotation per variant is used.
`columns_keep_snv`  | Columns to keep in the `SNV`/`TP53` sheets, selected from `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `FILTER`, `CALLER`, the extracted `vep_info_fields` (`SYMBOL` is renamed to `GENE`), and the extracted `format_fields`.
`snvs_remove`       | `Consequence` values to drop from the SNV sheets (e.g. `synonymous_variant`).
`sv_min_len`        | Minimum absolute `SVLEN` (in bp) for a chr14 SV to be included in the `SV` sheet. Defaults to `1000`.
`columns_drop_tn`   | Columns to drop from the `Translocations` sheet.
`columns_drop_sv`   | Columns to drop from the `SV` sheet.
`sv_col_max_widths` | Optional per-column width caps (in characters) for the `Translocations`/`SV` sheets, e.g. `{ALT: 30}`.

Example:

```yaml
format_fields:
  - GT
  - AF
  - DP
vep_info_fields:
  - Consequence
  - SYMBOL
  - HGVSc
  - HGVSp
columns_keep_snv:
  - CHROM
  - POS
  - REF
  - ALT
  - SYMBOL
  - Consequence
  - HGVSc
  - HGVSp
  - GT
  - AF
  - DP
snvs_remove:
  - synonymous_variant
sv_min_len: 1000
columns_drop_tn:
  - COVERAGE
columns_drop_sv:
  - COVERAGE
sv_col_max_widths:
  ALT: 30
```

### Gene list

`genes` under [`compile_xlsx_report`](/softwares/#compile_xlsx_report) can point to a BED file. If set, the `SNV` sheet
is restricted to variants whose gene (4th column of the BED file) is present in the file. This does not affect the
`TP53` sheet.

### Software versions

`software_versions` is passed as a dict of tool name to container image string (e.g.
`docker://hydragenetics/severus:1.5`), and the version is parsed from the tag after the last `:`. This is populated
automatically by the rule from other modules' `container` config values, and is not something that needs to be set
directly under `compile_xlsx_report`.
