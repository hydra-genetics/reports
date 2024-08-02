#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import yaml
import json
import csv


def metadata(sample, analysis_date, pipeline):
    # Uri Still in progress
    analysis_date = analysis_date
    pipeline_name = list(pipeline.keys())[0]
    pipeline = dict(name=pipeline_name,
                    version=pipeline[pipeline_name]['version'],
                    uri="")

    if tc != []:
        tc_pathology = tc[0]
        tc_purecn = tc[1]
    else:
        tc_pathology = 'NA'
        tc_purecn = 'NA'

    dict_ = dict(sample=sample,
                 analysis_date=analysis_date,
                 pipeline=pipeline,
                 tc_pathology=tc_pathology,
                 tc_purecn=tc_purecn)
    return dict_


def table(file, name, type, description, nav_header):
    with open(file) as f:
        list_ = []
        csv_table = csv.reader(f, delimiter='\t')
        header = next(csv_table)
        table = zip(*csv_table)
        for tupl in table:
            for val in tupl:
                list_.append(val)
        table_dict = dict(zip(header, list_))
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=[table_dict],
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
    if files['nav_header']:
        nav_header = files['nav_header']
    else:
        nav_header = ''
    return nav_header


def tmb(file, name, type, description, nav_header):
    with open(file) as f:
        tmb_report = f.readlines()
        tmb_value = tmb_report[0]
        nr_of_variants = tmb_report[1]
        tmb_results = tmb_value + '<br>' + nr_of_variants
    dict_ = dict(name=name,
                 type=type,
                 description=description,
                 value=tmb_results,
                 nav_header=nav_header)
    return dict_


def generate_json(output_files, sample, analysis_date, pipeline):
    with open(output_files, 'r') as f:
        output_files = yaml.safe_load(f)
    results = []
    for d in output_files['files']:
        path = d['input']
        sample_path = path.replace('{sample}_{type}', sample)
        nav_header = check_nav_header(d)
        if 'TMB' in d['input']:
            results1 = tmb(sample_path, d['name'], d['type'], d['description'],
                           nav_header)
        else:
            if d['type'] == 'table':
                results1 = table(sample_path, d['name'], d['type'],
                                 d['description'], nav_header)
            if d['type'] == 'file_table':
                results1 = file_table(sample_path, d['name'], d['type'],
                                      d['description'], nav_header)
            if d['type'] == 'file_link':
                results1 = file_table(sample_path, d['name'], d['type'],
                                      d['description'], nav_header)
            if d['type'] == 'multiqc':
                results1 = multiqc(sample_path, d['name'], d['type'],
                                   d['description'], ["table"], nav_header)
            if d["type"] == "image":
                results1 = file_table(sample_path, d['name'], d['type'],
                                      d['description'], nav_header)
        results.append(results1)
    meta_data = metadata(sample, analysis_date, pipeline)
    json_file = dict(meta_data, results=results)
    return json_file


if __name__ == "__main__":
    json_filename = snakemake.input.files
    output = snakemake.output.json
    output_files = snakemake.input.output_files
    tc = snakemake.params.tc
    sample = snakemake.params.sample
    analysis_date = snakemake.params.analysis_date
    pipeline = snakemake.params.pipeline_version
    json_file = generate_json(output_files, sample, analysis_date, pipeline)
    with open(output, 'w') as f:
        print(json.dumps(json_file), file=f)
