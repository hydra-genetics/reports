#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import yaml
import json
import datetime


def metadata(sample, pipeline, tc, units, reference_genome):
    # Uri Still in progress
    analysis_date = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    pipeline_name = list(pipeline.keys())[0]
    pipeline = dict(name=pipeline_name,
                    version=pipeline[pipeline_name]["version"],
                    reference_genome=reference_genome,
                    uri="")

    tc_pathology = tc[0]
    tc_purecn = tc[1]

    units = extract_sample_from_units(units, sample)
    dict_ = dict(sample=sample,
                 analysis_date=analysis_date,
                 pipeline=pipeline,
                 tc_pathology=tc_pathology,
                 tc_purecn=tc_purecn,
                 units=units)
    return dict_


def extract_sample_from_units(units, sample):
    sample = sample.split('_')[0]
    units = units.reset_index(drop=True)
    units = units[units['sample'] == sample]
    units_json = units.to_dict(orient='records')
    return units_json


def plain_text(file, name, type, description, nav_header):
    with open(file, "r") as f:
        plain_text = f.readlines()

    new_lines = []
    for val in plain_text:
        if "\n" in val:
            val = val.replace("\n", "")
        new_lines.append(val)
    plain_text = "<br>".join(new_lines)
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=plain_text,
                 nav_header=nav_header)
    return dict_


def table(file, name, type, description, nav_header):
    with open(file, "r") as f:
        json_file = json.load(f)

    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=json_file,
                 nav_header=nav_header)
    return dict_


def file_table(file, name, type, description, nav_header):
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=file,
                 nav_header=nav_header)
    return dict_


def multiqc(file, name, type, description, sections, nav_header):
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=file,
                 sections=sections,
                 nav_header=nav_header)
    return dict_


def check_nav_header(files):
    if files["nav_header"]:
        nav_header = files["nav_header"]
    else:
        nav_header = ""
    return nav_header


def tmb(file, name, type, description, nav_header):
    with open(file) as f:
        tmb_report = f.readlines()
        tmb_value = tmb_report[0]
        nr_of_variants = tmb_report[1]
        tmb_results = tmb_value + "<br>" + nr_of_variants
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=tmb_results,
                 nav_header=nav_header)
    return dict_


def generate_json(output_files, sample, pipeline, tc, units, reference_genome):
    with open(output_files, "r") as f:
        output_files = yaml.safe_load(f)
    results = []
    for d in output_files["files"]:
        path = d["input"]
        sample_path = path.replace("{sample}_{type}", sample)
        nav_header = check_nav_header(d)
        if "TMB" in d["input"]:
            results1 = tmb(sample_path, d["name"], d["type"], d["description"],
                           nav_header)
        else:
            if d["type"] == "table":
                results1 = table(sample_path, d["name"], d["type"],
                                 d["description"], nav_header)
            if d["type"] == "file_table":
                results1 = file_table(sample_path, d["name"], d["type"],
                                      d["description"], nav_header)
            if d["type"] == "large_file_table":
                results1 = file_table(sample_path, d["name"], d["type"],
                                      d["description"], nav_header)
            if d["type"] == "file_link":
                results1 = file_table(sample_path, d["name"], d["type"],
                                      d["description"], nav_header)
            if d["type"] == "multiqc":
                results1 = multiqc(sample_path, d["name"], d["type"],
                                   d["description"], ["table"], nav_header)
            if d["type"] == "image":
                results1 = file_table(sample_path, d["name"], d["type"],
                                      d["description"], nav_header)
            if d["type"] == "plain_text":
                results1 = plain_text(sample_path, d["name"], d["type"],
                                      d["description"], nav_header)
        results.append(results1)
    meta_data = metadata(sample, pipeline, tc, units, reference_genome)
    json_file = dict(meta_data, results=results)
    return json_file


if __name__ == "__main__":
    json_filename = snakemake.input.files
    output = snakemake.output.json
    output_files = snakemake.input.output_files
    tc = snakemake.params.tc
    sample = snakemake.params.sample
    pipeline = snakemake.params.pipeline_version
    units = snakemake.params.units
    reference_genome = snakemake.params.reference_genome
    json_file = generate_json(output_files, sample, pipeline, tc, units,
                              reference_genome)

    with open(output, "w") as f:
        print(json.dumps(json_file), file=f)
