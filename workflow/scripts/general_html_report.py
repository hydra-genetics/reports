from jinja2 import Template
import json
from jsonschema import validate
import time
import yaml


def validate_dict(d: dict, schema_path: str):
    with open(schema_path) as f:
        validate(instance=d, schema=yaml.safe_load(f))


def generate_report(template_filename: str, config: dict, final_directory_depth: int):
    with open(template_filename) as f:
        template = Template(source=f.read())

    if final_directory_depth != 0:
        for d in config["file_links"]:
            d["uri"] = final_directory_depth * "../" + d["uri"]

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
