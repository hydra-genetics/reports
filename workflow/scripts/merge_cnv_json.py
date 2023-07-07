from collections import defaultdict
from collections.abc import Generator
from dataclasses import dataclass
import json
from pathlib import Path
import pysam
import sys
from typing import Union
import cyvcf2


@dataclass
class CNV:
    caller: str
    chromosome: str
    genes: list
    start: int
    length: int
    type: str
    copy_number: float
    baf: float

    def end(self):
        return self.start + self.length - 1

    def overlaps(self, other):
        return self.chromosome == other.chromosome and (
            # overlaps in the beginning, or self contained in other
            (self.start >= other.start and self.start <= other.end())
            or
            # overlaps at the end, or self contained in other
            (self.end() >= other.start and self.end() <= other.end())
            or
            # other is contained in self
            (other.start >= self.start and other.end() <= self.end())
        )

    def __hash__(self):
        return hash(f"{self.caller}_{self.chromosome}:{self.start}-{self.end()}_{self.copy_number}")


def parse_fai(filename, skip=None):
    with open(filename) as f:
        for line in f:
            chrom, length = line.strip().split()[:2]
            if skip is not None and chrom in skip:
                continue
            yield chrom, int(length)


def parse_annotation_bed(filename, skip=None):
    with open(filename) as f:
        for line in f:
            chrom, start, end, name = line.strip().split()[:4]
            if skip is not None and chrom in skip:
                continue
            yield chrom, int(start), int(end), name


def get_vaf(vcf_filename, skip=None):
    vcf = cyvcf2.VCF(vcf_filename)
    for variant in vcf:
        if skip is not None and variant.CHROM in skip:
            continue
        yield variant.CHROM, variant.POS, variant.INFO.get("AF", None)


def get_cnvs(vcf_filename, skip=None):
    cnvs = defaultdict(lambda: defaultdict(list))
    vcf = pysam.VariantFile(vcf_filename)
    for variant in vcf.fetch():
        if skip is not None and variant.chrom in skip:
            continue
        caller = variant.info.get("CALLER")
        if caller is None:
            raise KeyError("could not find caller information for variant, has the vcf been annotated?")
        genes = variant.info.get("Genes")
        if genes is None:
            continue
        if isinstance(genes, str):
            genes = [genes]
        cnv = CNV(
            caller,
            variant.chrom,
            sorted(genes),
            variant.pos,
            variant.info.get("SVLEN"),
            variant.info.get("SVTYPE"),
            variant.info.get("CORR_CN"),
            variant.info.get("BAF"),
        )
        cnvs[variant.chrom][caller].append(cnv)
    return cnvs


def merge_cnv_dicts(dicts, vaf, annotations, chromosomes, filtered_cnvs, unfiltered_cnvs):
    callers = list(map(lambda x: x["caller"], dicts))
    caller_labels = dict(
        cnvkit="cnvkit",
        gatk="GATK",
    )
    cnvs = {}
    for chrom, chrom_length in chromosomes:
        cnvs[chrom] = dict(
            chromosome=chrom,
            label=chrom,
            length=chrom_length,
            vaf=[],
            annotations=[],
            callers={c: dict(name=c, label=caller_labels.get(c, c), ratios=[], segments=[], cnvs=[]) for c in callers},
        )

    for a in annotations:
        for item in a:
            cnvs[item[0]]["annotations"].append(
                dict(
                    start=item[1],
                    end=item[2],
                    name=item[3],
                )
            )

    if vaf is not None:
        for v in vaf:
            cnvs[v[0]]["vaf"].append(
                dict(
                    pos=v[1],
                    vaf=v[2],
                )
            )

    # Iterate over the unfiltered CNVs and pair them according to overlap.
    for uf_cnvs, f_cnvs in zip(unfiltered_cnvs, filtered_cnvs):
        for chrom, cnvdict in uf_cnvs.items():
            callers = sorted(list(cnvdict.keys()))
            first_caller = callers[0]
            rest_callers = callers[1:]

            # Keep track of added CNVs on a chromosome to avoid duplicates
            added_cnvs = set()

            for cnv1 in cnvdict[first_caller]:
                pass_filter = False

                if cnv1 in f_cnvs[chrom][first_caller]:
                    # The CNV is part of the filtered set, so all overlapping
                    # CNVs should pass the filter.
                    pass_filter = True

                cnv_group = [cnv1]
                for caller2 in rest_callers:
                    for cnv2 in cnvdict[caller2]:
                        if cnv1.overlaps(cnv2):
                            # Add overlapping CNVs from other callers
                            cnv_group.append(cnv2)

                            if cnv2 in f_cnvs[chrom][caller2]:
                                # If the overlapping CNV is part of the filtered
                                # set, the whole group should pass the filter.
                                pass_filter = True

                for c in cnv_group:
                    if c in added_cnvs:
                        continue
                    cnvs[c.chromosome]["callers"][c.caller]["cnvs"].append(
                        dict(
                            genes=c.genes,
                            start=c.start,
                            length=c.length,
                            type=c.type,
                            cn=c.copy_number,
                            baf=c.baf,
                            passed_filter=pass_filter,
                        )
                    )
                    added_cnvs.add(c)

    for d in dicts:
        for r in d["ratios"]:
            cnvs[r["chromosome"]]["callers"][d["caller"]]["ratios"].append(
                dict(
                    start=r["start"],
                    end=r["end"],
                    log2=r["log2"],
                )
            )
        for s in d["segments"]:
            cnvs[s["chromosome"]]["callers"][d["caller"]]["segments"].append(
                dict(
                    start=s["start"],
                    end=s["end"],
                    log2=s["log2"],
                )
            )

    for v in cnvs.values():
        v["callers"] = list(v["callers"].values())

    return list(cnvs.values())


def main():
    log = Path(snakemake.log[0])

    logfile = open(log, "w")
    sys.stdout = sys.stderr = logfile

    annotation_beds = snakemake.input["annotation_bed"]
    fasta_index_file = snakemake.input["fai"]
    germline_vcf = snakemake.input["germline_vcf"]
    json_files = snakemake.input["json"]
    filtered_cnv_vcf_files = snakemake.input["filtered_cnv_vcfs"]
    cnv_vcf_files = snakemake.input["cnv_vcfs"]

    if len(germline_vcf) == 0:
        germline_vcf = None

    output_file = snakemake.output["json"]

    skip_chromosomes = snakemake.params["skip_chromosomes"]

    cnv_dicts = []
    for fname in json_files:
        with open(fname) as f:
            cnv_dicts.append(json.load(f))

    fai = parse_fai(fasta_index_file, skip_chromosomes)
    vaf = None
    if germline_vcf is not None:
        vaf = get_vaf(germline_vcf, skip_chromosomes)
    annotations = []
    for filename in annotation_beds:
        annotations.append(parse_annotation_bed(filename, skip_chromosomes))

    filtered_cnv_vcfs = []
    unfiltered_cnv_vcfs = []
    for f_vcf, uf_vcf in zip(filtered_cnv_vcf_files, cnv_vcf_files):
        filtered_cnv_vcfs.append(get_cnvs(f_vcf, skip_chromosomes))
        unfiltered_cnv_vcfs.append(get_cnvs(uf_vcf, skip_chromosomes))

    cnvs = merge_cnv_dicts(cnv_dicts, vaf, annotations, fai, filtered_cnv_vcfs, unfiltered_cnv_vcfs)

    with open(output_file, "w") as f:
        print(json.dumps(cnvs), file=f)


if __name__ == "__main__":
    main()
