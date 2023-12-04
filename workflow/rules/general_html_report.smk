__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"


rule general_html_report:
    input:
        config_schema=workflow.source_path("../schemas/general_html_report_json.schema.yaml"),
        html_template=workflow.source_path("../templates/general_html_report/index.html"),
        json=config.get("general_html_report", {}).get("json"),
    output:
        html="reports/general_html_report/{sample}_{type}.general_report.html",
    params:
        final_directory_depth=config.get("general_html_report", {}).get("final_directory_depth", 1),
        extra=config.get("general_html_report", {}).get("extra", ""),
    log:
        "reports/general_html_report/{sample}_{type}.general_report.log",
    benchmark:
        repeat(
            "reports/general_html_report/{sample}_{type}.output.benchmark.tsv",
            config.get("general_html_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("general_html_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("general_html_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("general_html_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("general_html_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("general_html_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("general_html_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("general_html_report", {}).get("container", config["default_container"])
    message:
        "{rule}: generate general html report from json config {input.json}"
    script:
        "../scripts/general_html_report.py"
