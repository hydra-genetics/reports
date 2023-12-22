__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"


rule cnv_html_report:
    input:
        json="reports/cnv_html_report/{sample}_{type}.{tc_method}.merged.json",
        html_template=workflow.source_path("../templates/cnv_html_report/index.html"),
        js_files=[
            workflow.source_path("../templates/assets/js/d3.v7.min.js"),
            workflow.source_path("../templates/cnv_html_report/01-chromosome-plot.js"),
            workflow.source_path("../templates/cnv_html_report/02-genome-plot.js"),
            workflow.source_path("../templates/cnv_html_report/03-results-table.js"),
            workflow.source_path("../templates/cnv_html_report/04-window-summary.js"),
            workflow.source_path("../templates/cnv_html_report/05-main.js"),
        ],
        css_files=[
            workflow.source_path("../templates/cnv_html_report/style.css"),
        ],
        tc_file=get_tc_file,
    output:
        html=temp("reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html"),
    params:
        include_table=config.get("cnv_html_report", {}).get("show_table", True),
        tc=get_tc,
        tc_method=lambda wildcards: wildcards.tc_method,
        include_cytobands=config.get("cnv_html_report", {}).get("cytobands", False),
    log:
        "reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html.log",
    benchmark:
        repeat(
            "reports/cnv_html_report/{sample}_{type}.{tc_method}.cnv_report.html.benchmark.tsv",
            config.get("cnv_html_report", {}).get("benchmark_repeats", 1),
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
        json=temp("reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json"),
    params:
        skip_chromosomes=config.get("reference", {}).get("skip_chrs"),
    log:
        "reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json.log",
    benchmark:
        repeat(
            "reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json.benchmark.tsv",
            config.get("cnv_json", {}).get("benchmark_repeats", 1),
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
        annotation_bed=config.get("merge_cnv_json", {}).get("annotations", []),
        germline_vcf=get_germline_vcf,
        filtered_cnv_vcfs=get_filtered_cnv_vcf,
        cnv_vcfs=get_unfiltered_cnv_vcf,
        cytobands=config.get("merge_cnv_json", {}).get("cytobands", []),
    output:
        json=temp("reports/cnv_html_report/{sample}_{type}.{tc_method}.merged.json"),
    params:
        skip_chromosomes=config.get("reference", {}).get("skip_chrs", []),
        cytobands=config.get("cnv_html_report", {}).get("cytobands", False),
    log:
        "reports/cnv_html_report/{sample}_{type}.{tc_method}.merged.json.log",
    benchmark:
        repeat(
            "reports/cnv_html_report/{sample}_{type}.{tc_method}.merged.json.benchmark.tsv",
            config.get("merge_cnv_json", {}).get("benchmark_repeats", 1),
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
