from jinja2 import Template
from pathlib import Path
import time


def get_sample_name(filename):
    return Path(filename).name.split(".")[0]


def create_report(template_filename, json_filename, css_files, js_files, show_table, tc, tc_method):
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
    json_filename = snakemake.input.json
    template_dir = Path(snakemake.input.template_dir)
    html_filename = snakemake.output.html

    html_template = template_dir / "index.html"
    css_files = sorted(template_dir.glob("*.css"))
    js_files = sorted(template_dir.glob("*.js"))

    report = create_report(
        html_template,
        json_filename,
        css_files,
        js_files,
        snakemake.params.include_table,
        snakemake.params.tc,
        snakemake.params.tc_method,
    )

    with open(html_filename, "w") as f:
        f.write(report)


if __name__ == "__main__":
    main()
