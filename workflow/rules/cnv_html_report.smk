__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"


rule cnv_html_report:
    input:
        json="reporting/cnv_html_report/{sample}_{type}.{tc_method}.merged.json",
    output:
        report="reporting/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html",
    params:
        extra=config.get("cnv_html_report", {}).get("extra", ""),
    log:
        "reporting/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html.log",
    benchmark:
        repeat(
            "reporting/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html.benchmark.tsv",
            config.get("cnv_html_report", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cnv_html_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_html_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_html_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_html_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_html_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_html_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_html_report", {}).get("container", config["default_container"])
    message:
        "{rule}: Compile a CNV HTML report for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/cnv_html_report.py"


rule cnv_json:
    input:
        ratios=get_cnv_ratios,
        segments=get_cnv_segments,
    output:
        json="reporting/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json",
    params:
        extra=config.get("cnv_json", {}).get("extra", ""),
    log:
        "reporting/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json.log",
    benchmark:
        repeat(
            "reporting/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json.benchmark.tsv",
            config.get("cnv_json", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cnv_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_json", {}).get("container", config["default_container"])
    message:
        "{rule}: Create JSON representation for CNV results from {wildcards.caller} "
        "for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/cnv_json.py"


rule merge_cnv_json:
    input:
        json=get_json_for_merge_cnv_json,
        fai=config.get("reference", {}).get("fai", ""),
    output:
        json="reporting/cnv_html_report/{sample}_{type}.{tc_method}.merged.json",
    params:
        extra=config.get("merge_cnv_json", {}).get("extra", ""),
    log:
        "reporting/cnv_html_report/{sample}_{type}.{tc_method}.merged.json.log",
    benchmark:
        repeat(
            "reporting/cnv_html_report/{sample}_{type}.{tc_method}.merged.json.benchmark.tsv",
            config.get("merge_cnv_json", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("merge_cnv_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("merge_cnv_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge_cnv_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge_cnv_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge_cnv_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge_cnv_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge_cnv_json", {}).get("container", config["default_container"])
    message:
        "{rule}: Merge CNV JSON data for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/merge_cnv_json.py"
