from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
import pysam
import sys
from typing import Dict, Generator, List, Union


@dataclass
class CNV:
    caller: str
    chromosome: str
    genes: list
    start: int
    length: int
    type: str
    cn: float
    baf: float
    passed_filter: bool = False

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
        return hash(f"{self.caller}_{self.chromosome}:{self.start}-{self.end()}_{self.cn:.2f}")

    def __eq__(self, other):
        return hash(self) == hash(other)


def normalize_chrom(chrom: str) -> str:
    """Ensure chromosome name starts with 'chr' exactly once."""
    c = str(chrom)
    if c.startswith("chr"):
        c = c[3:]
    return f"chr{c}"


def parse_fai(filename, skip=None):
    with open(filename) as f:
        for line in f:
            chrom, length = line.strip().split()[:2]
            chrom = normalize_chrom(chrom)
            if skip is not None and chrom in skip:
                continue
            yield chrom, int(length)


def annotation_parser():
    parsed_annotations = set()

    def parse_annotation_bed(filename, skip=None):
        with open(filename) as f:
            for line in f:
                chrom, start, end, name = line.strip().split()[:4]
                chrom = normalize_chrom(chrom)
                if skip is not None and chrom in skip:
                    continue
                if (name, chrom, start, end) in parsed_annotations:
                    print(f"Warning: duplicate annotation {name} {chrom}:{start}-{end}", file=sys.stderr)
                    print(f"Warning: skipping {name} {chrom}:{start}-{end} in {filename}", file=sys.stderr)
                    continue
                parsed_annotations.add((name, chrom, start, end))
                yield chrom, int(start), int(end), name

    return parse_annotation_bed


def parse_cytobands(filename, cytoband_colors, cytoband_centromere="acen", skip=None):
    cytobands = defaultdict(list)
    with open(filename) as f:
        for line in f:
            chrom, start, end, name, giemsa = line.strip().split()
            chrom = normalize_chrom(chrom)
            if skip is not None and chrom in skip:
                continue
            cytobands[chrom].append(
                {
                    "name": name,
                    "start": int(start),
                    "end": int(end),
                    "direction": "none",
                    "giemsa": giemsa,
                    "color": cytoband_colors.get(giemsa, "#ff0000"),
                }
            )

    for k, v in cytobands.items():
        cytobands[k] = sorted(v, key=lambda x: x["start"])
        centromere_index = [i for i, x in enumerate(cytobands[k]) if x["giemsa"] == cytoband_centromere]

        if len(centromere_index) > 0 and len(centromere_index) != 2:
            print(
                f"error: chromosome {k} does not have 0 or 2 centromere bands, " f"found {len(centromere_index)}", file=sys.stderr
            )
            sys.exit(1)
        elif len(centromere_index) == 0:
            continue

        cytobands[k][centromere_index[0]]["direction"] = "right"
        cytobands[k][centromere_index[1]]["direction"] = "left"

    return cytobands


def parse_ref_genes(filename, skip=None):
    """Parse UCSC refGene.txt format to create a gene search index.
    
    Args:
        filename: Path to refGene.txt file
        skip: Set of chromosomes to skip
        
    Returns:
        Dict mapping gene names to {chrom, start, end}
    """
    if not filename:
        return {}
        
    genes = defaultdict(lambda: {"chrom": None, "start": float("inf"), "end": float("-inf")})
    
    with open(filename) as f:
        # Standard UCSC refGene.txt format (tab-separated):
        # With bin column: 0:bin, 1:name, 2:chrom, 4:txStart, 5:txEnd, 12:name2
        # Without bin column: 0:name, 1:chrom, 3:txStart, 4:txEnd, 11:name2
        
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 12:
                # Try fallback to any whitespace split if tabs aren't working as expected
                parts = line.strip().split()
                if len(parts) < 12:
                    continue
            
            # Detect layout
            if parts[0].isdigit() and len(parts) >= 13:
                # Likely has bin column
                chrom_idx = 2
                start_idx = 4
                end_idx = 5
                name_idx = 12
            elif parts[1].startswith("chr") or parts[1].isdigit() or parts[1].startswith(("X", "Y", "M")):
                # Likely no bin column
                chrom_idx = 1
                start_idx = 3
                end_idx = 4
                name_idx = 11
            else:
                # Fallback to defaults
                chrom_idx = 2
                start_idx = 4
                end_idx = 5
                name_idx = 12

            if name_idx >= len(parts):
                continue
                
            chrom = normalize_chrom(parts[chrom_idx])
            if skip is not None and chrom in skip:
                continue
                
            try:
                txStart = int(parts[start_idx])
                txEnd = int(parts[end_idx])
                name = parts[name_idx]
                
                if not name:
                    continue

                # Update gene extent (merge overlapping transcripts)
                if genes[name]["chrom"] is None:
                    genes[name]["chrom"] = chrom
                
                # Only merge if on same chromosome
                if genes[name]["chrom"] == chrom:
                    genes[name]["start"] = min(genes[name]["start"], txStart)
                    genes[name]["end"] = max(genes[name]["end"], txEnd)
            except (ValueError, IndexError):
                continue

    # Convert to standard dict and cleanup
    results = {k: v for k, v in genes.items() if v["chrom"] is not None}
    print(f"Parsed {len(results)} genes from {filename}", file=sys.stderr)
    return results


def get_vaf(vcf_filename: Union[str, bytes, Path], skip=None) -> Generator[tuple, None, None]:
    if skip is None:
        skip = []
    vcf = pysam.VariantFile(str(vcf_filename))
    for variant in vcf.fetch():
        chrom = normalize_chrom(variant.chrom)
        if chrom in skip:
            continue
        vaf = variant.info.get("AF")
        if isinstance(vaf, float):
            yield chrom, variant.pos, vaf
        elif vaf is not None:
            for f in vaf:
                yield chrom, variant.pos, f


def get_cnvs(vcf_filename, skip=None) -> Dict[str, Dict[str, List[CNV]]]:
    cnvs = defaultdict(lambda: defaultdict(list))
    vcf = pysam.VariantFile(vcf_filename)
    for variant in vcf.fetch():
        chrom = normalize_chrom(variant.chrom)
        if skip is not None and chrom in skip:
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
            chrom,
            sorted(genes),
            variant.pos,
            variant.info.get("SVLEN"),
            variant.info.get("SVTYPE"),
            variant.info.get("CORR_CN"),
            variant.info.get("BAF"),
        )
        cnvs[chrom][caller].append(cnv)
    return cnvs


def filter_chr_cnvs(unfiltered_cnvs: Dict[str, List[CNV]], filtered_cnvs: Dict[str, List[CNV]]) -> Dict[str, List[Dict]]:
    if len(unfiltered_cnvs) == 0:
        return {}
    callers = sorted(list(unfiltered_cnvs.keys()))
    first_caller = callers[0]
    rest_callers = callers[1:]

    # Keep track of added CNVs (and filter status) on a chromosome to avoid duplicates
    added_cnvs = {}

    for cnv1 in unfiltered_cnvs[first_caller]:
        pass_filter = False

        if cnv1 in filtered_cnvs.get(first_caller, []):
            # The CNV is part of the filtered set, so all overlapping
            # CNVs should pass the filter.
            pass_filter = True

        cnv_group = [cnv1]
        for caller2 in rest_callers:
            for cnv2 in unfiltered_cnvs[caller2]:
                if cnv1.overlaps(cnv2):
                    # Add overlapping CNVs from other callers
                    cnv_group.append(cnv2)

                    if cnv2 in filtered_cnvs.get(caller2, []):
                        # If the overlapping CNV is part of the filtered
                        # set, the whole group should pass the filter.
                        pass_filter = True

        for c in cnv_group:
            # Track CNV filter status.
            # If a CNV that was previously set to not pass the filter,
            # but later passes due to another comparison, then update
            # the filter status.
            if c not in added_cnvs or pass_filter:
                added_cnvs[c] = pass_filter

    cnvs = defaultdict(list)
    for c, pass_filter in added_cnvs.items():
        c.passed_filter = pass_filter
        cnvs[c.caller].append(c)
    return cnvs


def sort_cnvs(cnvs: List[CNV]) -> List[CNV]:
    def cnv_key_func(x):
        chrom_id = x.chromosome
        if x.chromosome.startswith("chr"):
            chrom_id = x.chromosome[3:]
        if chrom_id.isdigit():
            chrom_id = f"{int(chrom_id):0>5}"
        return [chrom_id, x.start]

    return sorted(cnvs, key=cnv_key_func)


def merge_cnv_calls(unfiltered_cnvs, filtered_cnvs):
    cnvs = []
    # Iterate over the filtered and unfiltered CNVs and pair them according to overlap.
    for uf_cnvs, f_cnvs in zip(unfiltered_cnvs, filtered_cnvs):
        for chrom in uf_cnvs.keys():
            merged_cnvs = filter_chr_cnvs(uf_cnvs[chrom], f_cnvs[chrom])
            for caller, m_cnvs in merged_cnvs.items():
                for c in m_cnvs:
                    if c not in cnvs:
                        cnvs.append(c)
                    else:
                        # If it already exists, make sure that it represents
                        # all genes and that the filtering status is updated
                        added_cnv = cnvs.pop(cnvs.index(c))
                        added_cnv.genes = list(set(added_cnv.genes + c.genes))
                        added_cnv.passed_filter = added_cnv.passed_filter or c.passed_filter
                        cnvs.append(added_cnv)
    return sort_cnvs(cnvs)


def merge_cnv_dicts(dicts, vaf, annotations, cytobands, chromosomes, filtered_cnvs, unfiltered_cnvs, gene_index=None):
    callers = list(map(lambda x: x["caller"], dicts))
    caller_labels = dict(
        cnvkit="cnvkit",
        gatk="GATK",
        jumble="jumble",
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

    for c in cytobands:
        cnvs[c]["cytobands"] = cytobands[c]

    if vaf is not None:
        for v in sorted(vaf, key=lambda x: x[1]):
            cnvs[v[0]]["vaf"].append(
                dict(
                    pos=v[1],
                    vaf=v[2],
                )
            )

    for cnv in merge_cnv_calls(unfiltered_cnvs, filtered_cnvs):
        cnvs[cnv.chromosome]["callers"][cnv.caller]["cnvs"].append(cnv)

    for d in dicts:
        for r in sorted(d["ratios"], key=lambda x: x["start"]):
            chrom = normalize_chrom(r["chromosome"])
            if chrom not in cnvs:
                continue
            cnvs[chrom]["callers"][d["caller"]]["ratios"].append(
                dict(
                    start=r["start"],
                    end=r["end"],
                    log2=r["log2"],
                )
            )
        for s in sorted(d["segments"], key=lambda x: x["start"]):
            chrom = normalize_chrom(s["chromosome"])
            if chrom not in cnvs:
                continue
            cnvs[chrom]["callers"][d["caller"]]["segments"].append(
                dict(
                    start=s["start"],
                    end=s["end"],
                    log2=s["log2"],
                )
            )

    for v in cnvs.values():
        v["callers"] = list(v["callers"].values())

    # Attach gene search index to first chromosome if available
    if gene_index and chromosomes:
        cnvs[chromosomes[0][0]]["gene_search_index"] = gene_index

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
    cytoband_file = snakemake.input["cytobands"]
    ref_genes_file = snakemake.input.get("ref_genes", "")

    if len(germline_vcf) == 0:
        germline_vcf = None

    output_file = snakemake.output["json"]

    skip_chromosomes = snakemake.params["skip_chromosomes"]
    if skip_chromosomes:
        skip_chromosomes = [normalize_chrom(c) for c in skip_chromosomes]
        
    show_cytobands = snakemake.params["cytobands"]

    cytoband_config = snakemake.config.get("merge_cnv_json", {}).get("cytoband_config", {}).get("colors", {})
    cytoband_centromere = "acen"
    cytoband_colors = {
        "gneg": cytoband_config.get("gneg", "#e3e3e3"),
        "gpos25": cytoband_config.get("gpos25", "#555555"),
        "gpos50": cytoband_config.get("gpos50", "#393939"),
        "gpos75": cytoband_config.get("gpos75", "#8e8e8e"),
        "gpos100": cytoband_config.get("gpos100", "#000000"),
        "acen": cytoband_config.get("acen", "#963232"),
        "gvar": cytoband_config.get("gvar", "#000000"),
        "stalk": cytoband_config.get("stalk", "#7f7f7f"),
    }

    cnv_dicts = []
    for fname in json_files:
        with open(fname) as f:
            cnv_dicts.append(json.load(f))

    fai = list(parse_fai(fasta_index_file, skip_chromosomes))
    vaf = None
    if germline_vcf is not None:
        vaf = get_vaf(germline_vcf, skip_chromosomes)

    annot_parser = annotation_parser()
    annotations = []
    for filename in annotation_beds:
        annotations.append(annot_parser(filename, skip_chromosomes))

    cytobands = []
    if cytoband_file and show_cytobands:
        cytobands = parse_cytobands(cytoband_file, cytoband_colors, cytoband_centromere, skip_chromosomes)

    # Parse gene index for comprehensive search
    gene_index = {}
    if ref_genes_file:
        gene_index = parse_ref_genes(ref_genes_file, skip_chromosomes)

    if len(filtered_cnv_vcf_files) != len(cnv_vcf_files):
        print(
            f"error: the number of unfiltered vcf files ({len(filtered_cnv_vcf_files)}) "
            f"must match the number of filtered vcf files ({len(cnv_vcf_files)})"
        )
        sys.exit(1)

    filtered_cnv_vcfs = []
    unfiltered_cnv_vcfs = []
    for f_vcf, uf_vcf in zip(filtered_cnv_vcf_files, cnv_vcf_files):
        filtered_cnv_vcfs.append(get_cnvs(f_vcf, skip_chromosomes))
        unfiltered_cnv_vcfs.append(get_cnvs(uf_vcf, skip_chromosomes))

    cnvs = merge_cnv_dicts(cnv_dicts, vaf, annotations, cytobands, fai, filtered_cnv_vcfs, unfiltered_cnv_vcfs, gene_index)

    with open(output_file, "w") as f:
        print(json.dumps(cnvs, default=vars), file=f)


if __name__ == "__main__":
    main()
