__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2026, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule reports_compile_xlsx_report:
    input:
        vcf_snv="snv_indels/whatshap_phase/{sample}_{type}.phased.include.panel.vep_annotated.vcf.gz",
        vcf_sv="cnv_sv/svdb_query/{sample}_{type}.pathology.svdb_query.vcf",
        vcf_cnv="cnv_sv/cnvkit_vcf/{sample}_{type}.pathology.annotate_cnv.germline.filter.cnv_filter.vcf",
    output:
        xlsx=temp("reports/xlsx_reports/{sample}_{type}_combined_report.xlsx"),
    params:
        filter_config=config.get("compile_xlsx_report", {}).get("filters", ""),
        genes_bed=config.get("compile_xlsx_report", {}).get("genes", ""),
        software_versions={
            "pbmm2": config.get("pbmm2_align", {}).get("container", ""),
            "vacmap": config.get("vacmap_align", {}).get("container", ""),
            "pbsv": config.get("pbsv_call", {}).get("container", ""),
            "sniffles2": config.get("sniffles2_call", {}).get("container", ""),
            "severus": config.get("severus_t_only", {}).get("container", ""),
            "cnvkit": config.get("cnvkit_batch", {}).get("container", ""),
            **(
                {
                    "clairs-to": config.get("clairs_to_call", {}).get("container", ""),
                    "deepsomatic": config.get("deepsomatic_t_only", {}).get("container", ""),
                }
                if config.get("use_deepsomatic", False)
                else {"clairs-to": config.get("clairs_to_call", {}).get("container", "")}
            ),
        },
    log:
        "reports/xlsx_reports/{sample}_{type}_combined_report.xlsx.log",
    benchmark:
        repeat(
            "reports/xlsx_reports/{sample}_{type}_combined_report.xlsx.benchmark.tsv",
            config.get("compile_xlsx_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("compile_xlsx_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("compile_xlsx_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("compile_xlsx_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("compile_xlsx_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("compile_xlsx_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("compile_xlsx_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("compile_xlsx_report", {}).get("container", config["default_container"])
    message:
        "{rule}: compile SNV and SV VCFs into Excel report for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/compile_xlsx_report.py"
