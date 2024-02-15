__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pathlib
import pandas as pd
import re
from typing import List, Union
import yaml
from snakemake.io import Wildcards
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("7.8.3")

### Set and validate config file

if not workflow.overwrite_configfiles:
    "At least one config file must be passed using --configfile/--configfiles, by command line or a profile!"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
            ]
        )

        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "_copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", f["name"].strip().lower()).split()))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


if len(workflow.modules) == 0:
    # Only generate copy-rules if the workflow is executed directly.
    generate_copy_rules(output_spec)


def get_cnv_callers(tc_method):
    for tcm in config.get("svdb_merge", {}).get("tc_method", []):
        if tcm["name"] == tc_method:
            return tcm["cnv_caller"]
    raise ValueError(f"no cnv caller config available for tc_method {tc_method}")


def get_json_for_merge_cnv_json(wildcards):
    callers = get_cnv_callers(wildcards.tc_method)
    return ["reports/cnv_html_report/{sample}_{type}.{caller}.{tc_method}.json".format(caller=c, **wildcards) for c in callers]


def get_cnv_ratios(wildcards):
    if wildcards.caller == "cnvkit":
        return "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr"

    if wildcards.caller == "gatk":
        return "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv"

    raise NotImplementedError(f"not implemented for caller {wildcards.caller}")


def get_cnv_segments(wildcards):
    if wildcards.caller == "cnvkit":
        return "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns"

    if wildcards.caller == "gatk":
        return "cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg"

    raise NotImplementedError(f"not implemented for caller {wildcards.caller}")


def get_germline_vcf(wildcards: Wildcards) -> List[Union[str, Path]]:
    return config.get("merge_cnv_json", {}).get("germline_vcf", [])


def get_germline_vcf_tbi(wildcards: Wildcards) -> List[Union[str, Path]]:
    vcf = get_germline_vcf(wildcards)
    if isinstance(vcf, list):
        return map(lambda x: x + ".tbi", vcf)
    else:
        return vcf + ".tbi"


def get_filtered_cnv_vcf(wildcards: Wildcards) -> List[Union[str, Path]]:
    if not config.get("cnv_html_report", {}).get("show_table", True):
        return []

    return config.get("merge_cnv_json", {}).get("filtered_cnv_vcfs", [])


def get_unfiltered_cnv_vcf(wildcards: Wildcards) -> List[Union[str, Path]]:
    if not config.get("cnv_html_report", {}).get("show_table", True):
        return []

    return config.get("merge_cnv_json", {}).get("unfiltered_cnv_vcfs", [])


def get_tc(wildcards):
    if wildcards.tc_method == "pathology":
        try:
            return get_sample(samples, wildcards)["tumor_content"]
        except KeyError:
            return None

    tc_file = get_tc_file(wildcards)

    if not os.path.exists(tc_file):
        return None

    with open(tc_file) as f:
        return f.read().strip()


def get_tc_file(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        return config["samples"]
    else:
        return f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
