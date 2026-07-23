# CNV report

The CNV HTML report provides an interactive, offline-capable visualization of copy number variant data. It is rendered using an HTML Canvas-based engine and does not require internet access.

## Input files

Required files are:

- log<sub>2</sub> ratios from each CNV caller
- Segments from each CNV caller

Optional input files are:

- Filtered and unfiltered VCFs that each contain calls from all included callers for generating a table of CNV calls
- BED files with annotations that should be added to plots
- Germline VCF file for displaying BAF in plots
- Cytoband definition for displaying in the chromosome plot

## Output files

- `reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html`

## Configuration

There are a couple of things that can be customised using the config file.

### Results table

The CNV results table contains CNVs that have been called by the pipeline. In order for the table to be included in the final report, `show_table` under [`cnv_html_report`](/softwares/#configuration) has to be `true`. If this is the case, then `unfiltered_cnv_vcfs` has to be defined under [`merge_cnv_json`](/softwares/#configuration_2).

Alongside the static `CN` column (the caller's corrected copy number at report-build time), the table includes an **Adjusted CN** column. This is recomputed live in the browser from each call's own raw segment log₂ ratio, using the same calculation the plots use, so it updates immediately as the baseline-offset, TC, and "Absolute copy number" controls are changed, without needing a report rebuild. Unlike the plots' Log2 ratio view, this column always applies purity correction — by default using the sample's estimated TC, not an assumed 100% purity — so it reflects a real adjusted copy number even before "Simulate purity" is checked; checking it and moving the TC slider lets you explore a different purity assumption. If a call has no matching segment for its position, this column shows "NA".

The **Type** column updates the same way when `table_filter_config` is set: it shows the name of the first group (see below) whose criteria the call's live Adjusted CN currently satisfies, title-cased (e.g. `amplification` → "Amplification"), or "Copy neutral" if no group matches. This replaces the caller's own static classification (whose vocabulary isn't standardized across callers — cnvkit reports `DUP`/`DEL`/`COPY_NORMAL`, GATK reports `<COPY_GAIN>`/`<COPY_LOSS>`/`<COPY_NORMAL>`) with one consistent, live vocabulary. Without `table_filter_config`, the column falls back to the caller's static classification, unchanged.

#### Filter toggle

The "filtered" checkbox above the table shows only calls that currently satisfy `table_filter_config` (see below), instead of a fixed pass/fail decision computed once at report-build time — so a call that only looks like a real amplification or loss under a different baseline/TC hypothesis becomes visible (or stops qualifying) as you adjust those controls, the same way the plots and the Adjusted CN column do. If a call from one caller doesn't currently qualify but an overlapping call from another caller does, the row is still shown (mirroring how overlapping calls are grouped elsewhere in the pipeline). If `table_filter_config` isn't set, the toggle has no criteria to filter by and shows every call.

`table_filter_config` (optional, under [`merge_cnv_json`](/softwares/#configuration_2)) is the path to a YAML file with any number of named groups — not limited to a fixed set of names. Each group is a tree of `any_of`/`all_of`/`not` nodes with `{field, operator, value}` leaves:

```yaml
amplification:
  all_of:
    - {field: adjusted_cn, operator: ">=", value: 6.0}
deletion:
  all_of:
    - {field: adjusted_cn, operator: "<=", value: 1.4}
```

A call's Type is the name of the **first** group (in the order they're defined in the file) whose criteria it satisfies, title-cased (e.g. `amplification` → "Amplification") unless the group sets an explicit `label` (see below); if none match, it's "Copy neutral". The filter toggle shows a call if it matches any group — unless that group sets `filter: false`, in which case the group still affects Type but is excluded from the filter toggle's decision. This is useful for an intermediate classification that doesn't itself warrant hard-filtering, e.g. a "duplication" group for gains too small to count as a full amplification:

```yaml
amplification:
  all_of:
    - {field: adjusted_cn, operator: ">=", value: 6.0}
duplication:
  filter: false  # Type-only: shows as "Duplication", but never makes the row visible under "filtered"
  all_of:
    - {field: adjusted_cn, operator: ">", value: 2.5}
    - {field: adjusted_cn, operator: "<", value: 6.0}
deletion:
  all_of:
    - {field: adjusted_cn, operator: "<", value: 1.4}
loh:
  label: "LOH"  # overrides the auto-generated "Loh" — useful for acronyms/exact wording
  all_of:
    - {field: adjusted_cn, operator: ">=", value: 1.4}
    - {field: adjusted_cn, operator: "<=", value: 2.5}
    - not:
        all_of:
          - {field: baf, operator: ">", value: 0.3}
          - {field: baf, operator: "<", value: 0.7}
```

The `loh` group above illustrates a real subtlety: a call can have a completely neutral copy number and still be clinically significant — a skewed BAF at CN≈2 indicates copy-neutral loss of heterozygosity (one allele lost, compensated by duplication of the other), distinct from an actual copy-number loss. Splitting these into separate named groups (rather than folding LOH detection into a single "deletion"/"loss" group) keeps the Type column from mislabeling a genuinely neutral-CN call as a deletion just because its BAF happens to be skewed.

Note that this `loh` group needs no caller-specific carve-out for callers that never measure BAF (e.g. Jumble): a leaf whose field has no data for a given call evaluates to "unknown" rather than false, and that unknown propagates through `not`/`any_of`/`all_of` instead of being silently treated as false. Concretely, `not({baf-in-neutral-range})` for a call with no BAF data evaluates to `not(unknown)`, which stays unknown — not `true` — so it does not, by itself, qualify the call as LOH. A plain two-valued evaluator would get this backwards (`not(false) = true`), incorrectly turning "we don't know the BAF" into "confirmed skewed" for every caller with missing BAF, not just one.

`operator` is one of `=`, `!=`, `<`, `>`, `<=`, `>=`. `field` may be:

- `adjusted_cn` — the call's live copy number (same value shown in the Adjusted CN column)
- `baf`, `length`, `caller` — always available, read from each call's VCF record
- any name listed in `extra_info_fields` (also under `merge_cnv_json`) — additional VCF INFO fields to extract per call and make available under that same name, for criteria that depend on something pipeline-specific (e.g. an artifact-database frequency annotation). Nothing beyond `adjusted_cn`/`baf`/`length`/`caller` is hardcoded, so this works the same way regardless of what that field happens to be called in a given pipeline's VCFs:

```yaml
merge_cnv_json:
  extra_info_fields:
    - Twist_AF
  table_filter_config: config/table_filter.yaml
```
```yaml
# config/table_filter.yaml
amplification:
  all_of:
    - {field: adjusted_cn, operator: ">=", value: 6.0}
    - {field: Twist_AF, operator: "<=", value: 0.15}
```

A gene-annotation-existence check isn't needed in this config — every call reaching the table already required a gene match to be included at all.

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

### Gene coloring

Gene coloring can be enabled in the chromosome plot by providing a CSV file with gene roles and colors via `cancer_genes` under [`merge_cnv_json`](/softwares/#configuration_2).

**Automatic Highlighting**:
If a gene is listed in the `cancer_genes` CSV and is also present in the `ref_genes` index (but not in the `annotations` BED files), its coordinates will be automatically fetched from the reference index and it will be highlighted on the chromosome plot. This allows for highlighting a large set of genes without manually creating individual BED entries. Genes defined in `annotations` BED files take precedence over auto-fetched coordinates.

Required CSV columns:
- `Gene`: Standard gene name (e.g., TP53)
- `Role`: Gene role (e.g., Oncogene, TSG) — used as a label in the report legend
- `Color`: Hex color code for the gene (e.g., #ff0000)

Example `cancer_genes.csv`:
```csv
Gene,Role,Color
TP53,Dual role (OG and / or TSG),#ee82ee
IKZF1,Tumor suppressor gene (TSG),#0000ff
```

The report will then include a toggle to apply these colors to the gene annotations in the chromosome plot.

### Manual tumor content adjustment

The report includes a slider that allows the user to manually override the estimated tumor cell content (TC) directly in the browser. The slider is **disabled by default** and only becomes active after enabling the **"Simulate purity"** checkbox in the chromosome view controls. Once enabled, adjusting the slider recalculates and redraws the expected copy number lines in real time without requiring a re-run of the pipeline. This is useful when the purity estimate is uncertain or when exploring alternative TC scenarios.

Enabling **"Simulate purity"** also switches both plots to the [Copy number view](#log2-ratio-vs-copy-number-view) by default (it can be switched back to log₂ ratio manually), and makes the **"Absolute copy number"** checkbox available. Enabling it snaps each segment line's recomputed copy number to the nearest whole number, which can make it easier to check the TC estimate against expected integer copy states. This does not affect the individual scatter points, only the segment lines. Unlike the view-mode toggle, "Absolute copy number" is only ever available while "Simulate purity" is checked.

### Log2 ratio vs. Copy number view

Both the chromosome and genome-wide plots can be switched between two Y-axis views:

- **Log2 ratio** (default) — the raw, purity-independent log₂ ratio. If "Simulate purity" is also enabled, values are TC-adjusted before being expressed as a log₂ ratio.
- **Copy number** — a linear, absolute copy-number scale (0–5 by default) instead of log₂ ratio. This view works even without "Simulate purity" enabled, in which case it assumes 100% purity; enabling "Simulate purity" makes the displayed numbers reflect the actual TC-adjusted copy number.

**Copy number is the stable anchor in both views**: adjusting the baseline offset slider never
moves segments or points — a given absolute copy number always renders at the same position, in
either view. What changes instead is how the axes are **labelled**:

- In **Log2 ratio** view, the primary (left) axis's fixed positions are relabelled so the "0" label
  moves to wherever the current baseline sits — its underlying scale and the plotted data never
  move. The secondary (right) copy-number axis keeps its standard, fixed tick *positions* (1, 2, 4,
  8...), but its printed labels rescale the same multiplicative way as Copy Number view's primary
  axis, so the two views read consistently with each other under a shifted baseline.
- In **Copy number** view, the primary axis's fixed whole-number positions are instead relabelled
  multiplicatively (a fixed row showing "2" at baseline 0 relabels to "2 × 2^baseline" once the
  baseline changes) — ticks are generally no longer round integers once the baseline is nonzero.

Small triangle markers on the left and right edges of the plot point at wherever the current
baseline sits (the row currently labelled "2" / "0"), while a thicker, full-width reference line
always marks the fixed position of true absolute copy number 2 (diploid) — independent of the
baseline adjustment, so there's a stable anchor to judge an adjusted baseline against.

The secondary axis on the right-hand side of each plot always shows copy number regardless of the
active view (in Copy number view it mirrors the primary axis exactly), and the BAF row is labelled
on both the left and right for readability.

### Y-axis zoom

The "Y-axis zoom" controls (+ / − / Reset) narrow or widen the visible range of the ratio panel
(log2 ratio or copy number, depending on the active view) around its current center, on both the
chromosome and genome-wide plots. This is independent of "Zoom to data extent" and of the X-axis
(genomic position) zoom, and doesn't affect the BAF panel, whose Y-range is always fixed to
`[0, 1]`. "Reset" restores the default (unzoomed) range.

### Baseline and tumor content

The report includes a baseline offset slider, driven by the mechanism described above.

Once a baseline offset is active (entered manually), the report remembers which raw log₂ value it currently treats as exactly 2 copies. Because the mapping between a raw (TC-diluted) log₂ value and its absolute copy number depends on TC once "Simulate purity" is in effect, adjusting the TC slider afterwards — or toggling "Simulate purity" itself, which switches the TC actually applied between 1 and the real value — automatically re-solves the baseline offset so that same reference point stays pinned at exactly 2 copies, rather than silently drifting out of sync with the newly-adjusted TC.

### Gene Focus

The "Gene Focus" checkbox (chromosome view) displays data points with equal spacing along the x-axis, regardless of their genomic position. This allows for visualizing the sequential order of probes rather than their physical distance, which can be useful when probes are unevenly distributed — e.g. to inspect a densely-probed gene region without it being compressed by the surrounding sparser intergenic distance.

### Default-checked controls

Several checkboxes can be configured to be checked by default when the report is opened, via `cnv_html_report` in the config file:

```yaml
cnv_html_report:
    default_simulate_purity: true
    default_absolute_copy_number: true
    default_gene_focus: true
    default_show_all_datapoints: true
```

- `default_simulate_purity` — checks "Simulate purity" (and, per the behavior described above, switches the default view to Copy number).
- `default_absolute_copy_number` — checks "Absolute copy number". Only takes effect if `default_simulate_purity` is also `true`, since the checkbox requires purity simulation to be active; if set without it, the checkbox is forced back to unchecked when the report loads.
- `default_gene_focus` — checks "Gene Focus".
- `default_show_all_datapoints` — checks "Show all data points".

All four default to `false`.

### Wide Plots

To support wider plots that stack vertically (instead of the default responsive layout), you can configure `wide_plot_width` under [`cnv_html_report`](/softwares/#configuration). This accepts an integer value in pixels.

```yaml
cnv_html_report:
    wide_plot_width: 2000
```

### Interactive visualization features

The report includes the following interactive features:

| Feature | Description |
|---|---|
| Chromosome plot | Log₂-ratio and BAF scatter plot per chromosome with zoom and pan |
| Genome-wide plot | Overview of copy number across all chromosomes |
| Linear chromosome view | Alternative linear view for each chromosome with per-caller toggle |
| Gene search | Search box to quickly navigate the plot to a specific gene |
| Log2 ratio / Copy number toggle | Switches the Y-axis between log₂ ratio and a linear, absolute copy-number scale; see [Log2 ratio vs. Copy number view](#log2-ratio-vs-copy-number-view) |
| Manual TC adjustment | Slider to override estimated tumor content and update copy number lines in real time; requires **Simulate purity** to be enabled first |
| Absolute copy number | Snaps segment lines to whole-number copy number; requires **Simulate purity** to be enabled first |
| Gene Focus | Displays data points with equal spacing along the x-axis instead of by genomic position |
| Gene color toggle | Toggle to apply per-gene role colors to annotated genes in the plot |
| Caller toggle | Switch between callers (CNVkit, GATK, Jumble) in the chromosome and genome plots |
| Adjusted CN column | Live-recomputed copy-number column in the results table, following the baseline/TC/rounding controls; see [Results table](#results-table) |


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
        ploidy_file=reports.get_ploidy_file,
```

The only template file that is strictly required is the html file `html_template`. Both `css_files` and `js_files` can be left out if you so wish, but the functionality will be severly limited without any javascript.
