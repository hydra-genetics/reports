from jinja2 import Template
import json
from jsonschema import validate
import time
import yaml


def validate_dict(d: dict, schema_path: str):
    with open(schema_path) as f:
        validate(instance=d, schema=yaml.safe_load(f))


def generate_report(template_filename: str, config: dict):
    with open(template_filename) as f:
        template = Template(source=f.read())

    return template.render(
        dict(
            metadata=dict(
                analysis_date=config["analysis_date"],
                report_date=time.strftime("%Y-%m-%d %H:%M", time.localtime()),
                sample=config["sample"],
            ),
            pipeline=config["pipeline"],
            file_links=config["file_links"],
        )
    )


def main():
    html_template = snakemake.input.html_template
    json_file = snakemake.input.json

    with open(json_file) as f:
        config = json.load(f)

    validate_dict(config, snakemake.input.config_schema)

    report_content = generate_report(html_template, config)

    with open(snakemake.output.html, "w") as f:
        f.write(report_content)


if __name__ == "__main__":
    main()
