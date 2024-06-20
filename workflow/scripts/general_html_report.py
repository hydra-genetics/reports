from jinja2 import Template
import json
from jsonschema import validate
from jsonschema.exceptions import ValidationError
import re
import sys
import time
import yaml
import pandas as pd


def validate_dict(d: dict, schema_path: str):
    with open(schema_path, 'r') as f:
        try:
            validate(instance=d, schema=yaml.safe_load(f))
        except ValidationError as ve:
            print(f"error: failed to validate general report config:",
                  file=sys.stderr)
            print(ve, file=sys.stderr)
            sys.exit(1)


def validate_table_data(table: list) -> bool:
    if len(table) == 0:
        raise ValueError("empty table")

    n_cols = len(table[0])
    cols = table[0].keys()

    for i, row in enumerate(table, start=1):
        if len(row) != n_cols:
            raise ValueError(
                f"expected {n_cols} columns in row {i}, found {len(row)}")
        if set(row.keys()) != set(cols):
            raise ValueError(
                f"expected columns {', '.join(repr(x) for x in cols)} "
                f"in row {i}, found {', '.join(repr(x) for x in row.keys())}")
    return True


def fix_relative_uri(uri: str, depth: int) -> str:
    if re.match("^([^/]+:|/)", uri):
        return uri
    return depth * "../" + uri


def parse_multiqc(d: dict):
    with open(d["value"]) as f:
        multiqc_dict = json.loads(f.read())

    multiqc_res = {}

    for s in d["sections"]:
        if s == "table":
            multiqc_res[s] = {}
            multiqc_res[s]["data"] = {}
            multiqc_res[s]["header"] = {}

            for h_section, d_section in zip(
                    multiqc_dict["report_general_stats_headers"],
                    multiqc_dict["report_general_stats_data"]):
                for k, v in h_section.items():
                    multiqc_res[s]["header"][k] = v
                for sample, cols in d_section.items():
                    if sample not in multiqc_res[s]["data"]:
                        multiqc_res[s]["data"][sample] = {}
                    for k, v in h_section.items():
                        if k not in cols:
                            continue
                        multiqc_res[s]["data"][sample][k] = cols[k]
    return multiqc_res


def navigation_bar(config: dict):
    headers = []
    for d in config["results"]:
        headers.append(d["nav_header"])
    headers_no_dup = list(set(headers))
    return headers_no_dup


def merge_json(config: dict, extra_config: dict):
    for d in extra_config['results']:
        config['results'].append(d)
    return config


def generate_report(template_filename: str, config: dict,
                    final_directory_depth: int, css_files: list,
                    navigation_bar: list) -> str:
    with open(template_filename) as f:
        template = Template(source=f.read())

    for d in config["results"]:
        if d["type"] == "file_link":
            d["value"] = fix_relative_uri(d["value"], final_directory_depth)
        if d["type"] == "table":
            validate_table_data(d["value"])
        if d["type"] == "image":
            d["value"] = fix_relative_uri(d["value"], final_directory_depth)

        if d["type"] == "multiqc":
            d["data"] = parse_multiqc(d)

        if d["type"] == "file_table":
            data = pd.read_csv(d['value'], sep='\t')
            if len(data.columns) == 1:
                data = pd.read_csv(d['value'], sep=',')
            d["data"] = data.to_html(index=False).replace('border="1"', '')

    css_string = ""
    for css_filename in css_files:
        with open(css_filename) as f:
            css_string += f.read()

    navigation_bar.sort()
    nav_bar_html = ''

    # Rewrite in javascript
    for header in navigation_bar:
        nav_bar_html += f'\t\t\t\t<button class="tablinks" onclick="openNav(event, \'{header}\')">{header}</button>\n'

    return template.render(
        dict(metadata=dict(
            analysis_date=config["analysis_date"],
            report_date=time.strftime("%Y-%m-%d %H:%M", time.localtime()),
            sample=config["sample"],
        ),
             pipeline=config["pipeline"],
             results=config["results"],
             css=css_string,
             nav_bar=nav_bar_html,
             nav_header=navigation_bar))


def main():
    html_template = snakemake.input.html_template
    json_file = snakemake.input.json
    additional_json_file = snakemake.input.additional_json
    css = snakemake.input.css_files
    config_schema = snakemake.input.config_schema
    final_directory_depth = snakemake.params.final_directory_depth

    with open(json_file) as f:
        config = json.load(f)

    if len(additional_json_file) == 0:
        config = config
    else:
        with open(additional_json_file) as f:
            additional_config = json.load(f)
            config = merge_json(config, additional_config)

    nav_bar = navigation_bar(config)
    validate_dict(config, schema_path=config_schema)
    report_content = generate_report(html_template, config,
                                     final_directory_depth, css, nav_bar)

    with open(snakemake.output.html, "w") as f:
        f.write(report_content)


if __name__ == "__main__":
    main()
