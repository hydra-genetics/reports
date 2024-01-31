import csv
from jinja2 import Template
from pathlib import Path
import sys
import time


def get_sample_name(filename):
    return Path(filename).name.split(".")[0]


def parse_table(table_def):
    tsv_path = table_def["path"].format(**snakemake.wildcards)
    with open(tsv_path) as f:
        table_data = list(csv.DictReader(f, delimiter="\t"))

    return {
        "name": table_def.get("name", ""),
        "description": table_def.get("description", ""),
        "header": list(table_data[0].keys()),
        "data": table_data,
    }


def create_report(template_filename, json_filename, css_files, js_files, show_table, extra_tables, tc, tc_method):
    with open(template_filename) as f:
        template = Template(source=f.read())

    with open(json_filename) as f:
        json_string = f.read()

    css_string = ""
    for css_filename in css_files:
        with open(css_filename) as f:
            css_string += f.read()

    js_string = ""
    for js_filename in js_files:
        with open(js_filename) as f:
            js_string += f.read()

    return template.render(
        dict(
            json=json_string,
            css=css_string,
            js=js_string,
            extra_tables=extra_tables,
            metadata=dict(
                date=time.strftime("%Y-%m-%d %H:%M", time.localtime()),
                sample=get_sample_name(json_filename),
                show_table=show_table,
                tc=tc,
                tc_method=tc_method,
            ),
        )
    )


def main():
    log = Path(snakemake.log[0])

    logfile = open(log, "w")
    sys.stdout = sys.stderr = logfile

    json_filename = snakemake.input.json
    html_template = Path(snakemake.input.html_template)
    html_filename = snakemake.output.html
    extra_tables_def = snakemake.params.extra_tables

    js_files = []
    if "js_files" in snakemake.input.keys():
        js_files = snakemake.input.js_files

    css_files = []
    if "css_files" in snakemake.input.keys():
        css_files = snakemake.input.css_files

    extra_tables = [parse_table(t) for t in extra_tables_def]

    report = create_report(
        html_template,
        json_filename,
        css_files,
        js_files,
        snakemake.params.include_table,
        extra_tables,
        snakemake.params.tc,
        snakemake.params.tc_method,
    )

    with open(html_filename, "w") as f:
        f.write(report)


if __name__ == "__main__":
    main()
