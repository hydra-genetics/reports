# Icons

The use of icons in the report templates is enabled by embedding a minimal subset of [Bootstrap Icons](https://icons.getbootstrap.com/). This can be found in the CSS file located at [`workflow/templates/assets/css/icons.css`](https://github.com/hydra-genetics/reports/tree/develop/workflow/templates/assets/css/icons.css). The font itself is embedded as a base64 encoded string within this CSS file, making sure that no internet connection is required in order to use this feature. Simply include this CSS file in the Snakemake rule for the report in question and integrate it in the template:

```python
rule my_report:
    input:
        css: [
            workflow.source_path("../templates/assets/icons/icons.css"),
        ]
```

## Adding icons

The utility [`scripts/bootstrap_icons_subset.sh`](https://github.com/hydra-genetics/reports/tree/develop/scripts/bootstrap_icons_subset.sh) can be used to generate new subsets of icons. The only input required is a comma-separated string of of unicode code points that should be included in the resulting font. A base64-encoded string representing the font will be output on stdout, while the cmap table and file information will be output on stderr.

To embed the font in a CSS file, use the snipped below and replace `<font string>` with the output from the script.

```css
@font-face {
  font-display: block;
  font-family: "bootstrap-icons";
  src: url("data:font/woff2;charset=utf-8;base64,<font string>")
    format("woff2");
}
```
