from jinja2 import Template
import json
from jsonschema import validate
import re
import time
import yaml


def validate_dict(d: dict, schema_path: str):
    with open(schema_path) as f:
        validate(instance=d, schema=yaml.safe_load(f))


def validate_table_data(table: list) -> bool:
    if len(table) == 0:
        raise ValueError("empty table")

    n_cols = len(table[0])
    cols = table[0].keys()

    for i, row in enumerate(table, start=1):
        if len(row) != n_cols:
            raise ValueError(f"expected {n_cols} columns in row {i}, found {len(row)}")
        if set(row.keys()) != set(cols):
            raise ValueError(
                f"expected columns {', '.join(repr(x) for x in cols)} "
                f"in row {i}, found {', '.join(repr(x) for x in row.keys())}"
            )

    return True


def fix_relative_uri(uri: str, depth: int) -> str:
    if re.match("^([^/]+:|/)", uri):
        return uri

    return depth * "../" + uri


def generate_report(template_filename: str, config: dict, final_directory_depth: int) -> str:
    with open(template_filename) as f:
        template = Template(source=f.read())

    for d in config["file_links"]:
        d["uri"] = fix_relative_uri(d["uri"], final_directory_depth)

    for d in config["results"]:
        if d["type"] == "table":
            validate_table_data(d["value"])
        if d["type"] == "image":
            d["value"] = fix_relative_uri(d["value"], final_directory_depth)

    return template.render(
        dict(
            metadata=dict(
                analysis_date=config["analysis_date"],
                report_date=time.strftime("%Y-%m-%d %H:%M", time.localtime()),
                sample=config["sample"],
            ),
            pipeline=config["pipeline"],
            file_links=config["file_links"],
            results=config["results"],
        )
    )


def main():
    html_template = snakemake.input.html_template
    json_file = snakemake.input.json
    final_directory_depth = snakemake.params.final_directory_depth

    with open(json_file) as f:
        config = json.load(f)

    validate_dict(config, snakemake.input.config_schema)

    report_content = generate_report(html_template, config, final_directory_depth)

    with open(snakemake.output.html, "w") as f:
        f.write(report_content)


if __name__ == "__main__":
    main()
