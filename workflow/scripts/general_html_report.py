from jinja2 import Template


def generate_report(template_filename, json):
    with open(template_filename) as f:
        template = Template(source=f.read())

    return template.render(dict(sample=snakemake.wildcards.sample))


def main():
    html_template = snakemake.input.html_template
    json = snakemake.input.json

    report_content = generate_report(html_template, json)

    with open(snakemake.output.html, "w") as f:
        f.write(report_content)


if __name__ == "__main__":
    main()
