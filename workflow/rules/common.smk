__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"

import csv
import itertools
import linecache
import numpy as np
import pathlib
import pandas as pd
import re
import sys
from typing import List, Union
import yaml
from snakemake.iocontainers import Wildcards
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import datetime
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.software_versions import get_pipeline_version

min_version("9.0.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str)

if units.platform.iloc[0] in ["PACBIO", "ONT"]:
    units = units.set_index(["sample", "type", "processing_unit", "barcode"], drop=False).sort_index()
else:  # assume that the platform Illumina data with a lane and flowcell columns
    units = units.set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()

validate(units, schema="../schemas/units.schema.yaml")

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")

pipeline_name = ""
pipeline_version = get_pipeline_version(workflow, pipeline_name=pipeline_name)


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),
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


class _IdentityLineMap(dict):
    """Line mapping for generated code, whose compiled and source line numbers are the same."""

    def __contains__(self, lineno):
        return True

    def __missing__(self, lineno):
        return lineno


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
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb={mem_mb}, '
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

    source = "\n".join(rulestrings)
    source_name = "copy_result_files"

    # Snakemake 9 derives rule.run_func_src for the "code" rerun trigger by looking the compiled
    # filename up in workflow.linemaps and reading the lines back through linecache. Neither knows
    # about source we compile ourselves, so register it in both. The generated source is its own
    # original, so compiled and source line numbers coincide.
    linecache.cache[source_name] = (len(source), None, source.splitlines(True), source_name)
    workflow.linemaps[source_name] = _IdentityLineMap()

    exec(compile(source, source_name, "exec"), workflow.globals)


with open(config["general_report"]) as f:
    if f.name.endswith(".yaml"):
        general_report = yaml.safe_load(f)


def is_imported_as_module():
    """Whether this workflow is being parsed as a module of a consuming workflow.

    Snakemake 7 registered modules in ``workflow.modules``, which stayed empty when the
    workflow ran directly. In Snakemake 9 that attribute is empty in both cases, so the
    check has moved to the workflow modifier: the top-level workflow is parsed with the
    base modifier, which has no parent, while a module is parsed with a child modifier
    whose ``parent_modifier`` is the importing workflow's.
    """
    modifier = getattr(workflow, "modifier", None)
    return getattr(modifier, "parent_modifier", None) is not None


if not is_imported_as_module():
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

    if wildcards.caller == "jumble":
        return "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.cnr"

    raise NotImplementedError(f"not implemented for caller {wildcards.caller}")


def get_cnv_segments(wildcards):
    if wildcards.caller == "cnvkit":
        return "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns"

    if wildcards.caller == "gatk":
        return "cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg"

    if wildcards.caller == "jumble":
        return "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.cns"

    raise NotImplementedError(f"not implemented for caller {wildcards.caller}")


def get_germline_vcf(wildcards: Wildcards) -> List[Union[str, Path]]:
    return config.get("merge_cnv_json", {}).get("germline_vcf", [])


def get_unfiltered_cnv_vcf(wildcards: Wildcards) -> List[Union[str, Path]]:
    if not config.get("cnv_html_report", {}).get("show_table", True):
        return []

    return config.get("merge_cnv_json", {}).get("unfiltered_cnv_vcfs", [])


def get_cytobands(wildcards: Wildcards) -> List[Union[str, Path]]:
    return config.get("merge_cnv_json", {}).get("cytobands", [])


def get_ref_genes(wildcards: Wildcards) -> List[Union[str, Path]]:
    return config.get("merge_cnv_json", {}).get("ref_genes", [])


def get_cancer_genes(wildcards: Wildcards) -> List[Union[str, Path]]:
    res = config.get("merge_cnv_json", {}).get("cancer_genes", [])
    if isinstance(res, str) and not res:
        return []
    return res


if not config.get("merge_cnv_json", {}).get("cancer_genes"):
    print("WARNING: merge_cnv_json: cancer_genes not specified. Gene coloring will be disabled in the report.")


def get_tc(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology_purecn":
        tc = ""
        tc_file = f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if os.path.exists(tc_file):
            with open(tc_file) as f:
                tc = f.read().strip()
        if tc == "" or float(tc) < 0.35:
            return get_sample(samples, wildcards)["tumor_content"]
        else:
            return tc
    elif tc_method == "pathology":
        return get_sample(samples, wildcards)["tumor_content"]
    else:
        tc_file = f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if not os.path.exists(tc_file):
            return -1
        else:
            with open(tc_file) as f:
                tc = f.read().strip()
                if tc == "":
                    return "0.2"
                else:
                    return tc


def get_tc_general_report(wildcards):
    if get_sample(samples, wildcards)["tumor_content"]:
        tc_pathology = get_sample(samples, wildcards)["tumor_content"]
    else:
        tc_pathology = "NA"
    tc_purecn = "NA"
    tc_file = f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
    if os.path.exists(tc_file):
        with open(tc_file) as f:
            tc_purecn = f.read().strip()
    return [tc_pathology, tc_purecn]


def get_tc_file(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        return config["samples"]
    else:
        return f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"


def get_ploidy_file(wildcards):
    ploidy_file = f"cnv_sv/purecn/{wildcards.sample}_{wildcards.type}.csv"
    if os.path.exists(ploidy_file):
        return ploidy_file
    return config["samples"]


def get_ploidy(wildcards):
    ploidy_file = f"cnv_sv/purecn/{wildcards.sample}_{wildcards.type}.csv"
    if not os.path.exists(ploidy_file):
        return None
    with open(ploidy_file) as f:
        reader = csv.DictReader(f)
        try:
            row = next(reader)
            return float(row["Ploidy"])
        except (StopIteration, ValueError, KeyError, TypeError):
            return None
