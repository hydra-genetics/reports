from collections import defaultdict
from dataclasses import dataclass
import json
import math
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
    while c.startswith("chr"):
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

    print(f"DEBUG: Opening ref_genes file: {filename}", file=sys.stderr)

    with open(filename) as f:
        # Standard UCSC refGene.txt format (tab-separated):
        # With bin column: 0:bin, 1:name, 2:chrom, 4:txStart, 5:txEnd, 12:name2
        # Without bin column: 0:name, 1:chrom, 3:txStart, 4:txEnd, 11:name2

        line_count = 0
        for line in f:
            line_count += 1
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
                # Likely has bin column: 0:bin, 1:name, 2:chrom, 4:txStart, 5:txEnd, 12:name2
                chrom_idx = 2
                start_idx = 4
                end_idx = 5
                name_idx = 12
            elif parts[1].startswith("chr") or parts[1].isdigit() or parts[1].startswith(("X", "Y", "M")):
                # Likely no bin column: 0:name, 1:chrom, 3:txStart, 4:txEnd, 11:name2
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

            # Skip header row if present
            if parts[start_idx] == "txStart" or parts[chrom_idx] == "chrom":
                continue

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
            except Exception as e:
                # Log first few errors to help debugging, but don't crash the whole script
                if line_count < 20:
                    print(f"DEBUG: Error parsing line {line_count}: {e}. Row content: {parts[:6]}...", file=sys.stderr)
                continue

    # Convert to standard dict and cleanup
    results = {k: v for k, v in genes.items() if v["chrom"] is not None}
    print(f"Parsed {len(results)} genes from {filename} into search index", file=sys.stderr)
    return results


def get_baf(vcf_filename: Union[str, bytes, Path], skip=None) -> Generator[tuple, None, None]:
    if skip is None:
        skip = []
    vcf = pysam.VariantFile(str(vcf_filename))
    variants = []
    # Fetch all records once since fetch() might be expensive or stateful
    records = list(vcf.fetch())
    found_any = False
    for sample in vcf.header.samples:
        vaf_count = 0
        for record in records:
            # Try BAF, then VAF, then AF from INFO
            baf = None
            if "BAF" in record.info: baf = record.info["BAF"]
            elif "VAF" in record.info: baf = record.info["VAF"]
            elif "AF" in record.info: baf = record.info["AF"]

            if baf is None:
                # Check sample
                s = record.samples[sample]
                if "BAF" in s: baf = s["BAF"]
                elif "VAF" in s: baf = s["VAF"]
                elif "AF" in s: baf = s["AF"]
                elif "AD" in s and "DP" in s:
                    ad = s["AD"]
                    dp = s["DP"]
                    if ad and dp and dp > 0 and len(ad) >= 2:
                        baf = ad[1] / dp

            if baf is not None:
                chrom = normalize_chrom(record.chrom)
                if skip is not None and chrom in skip:
                    continue
                if isinstance(baf, (list, tuple)):
                    for bb in baf:
                        if bb is not None: variants.append((chrom, record.pos, float(bb)))
                else:
                    variants.append((chrom, record.pos, float(baf)))
                vaf_count += 1
        
        if vaf_count > 0:
            found_any = True
            break # Found data for a sample

    vcf.close()
    return variants


def bin_baf(baf_list, poi_regions, roi_flank_size_bp=10000, target_data_points=10000, roi_budget_fraction=0.5, roi_resolution_factor=10):
    if not baf_list:
        return []

    print(f"DEBUG: bin_baf - Input points: {len(baf_list)}, Target: {target_data_points}", file=sys.stderr)

    # Ensure baf_list is sorted by position
    baf_list.sort(key=lambda x: (x[0], x[1]))

    # Group ROI by chromosome and cover the entire region
    poi_by_chrom = defaultdict(list)
    for s in poi_regions:
        full_region = False
        if isinstance(s, dict):
            start, end, chrom = s["start"], s["end"], s["chromosome"]
            full_region = s.get("full_region", False)
        else:
            chrom, start, end = s[0], s[1], s[2]
            if len(s) > 4: full_region = s[4]

        if full_region:
            poi_by_chrom[chrom].append((start - roi_flank_size_bp, end + roi_flank_size_bp))
        else:
            poi_by_chrom[chrom].append((start - roi_flank_size_bp, start + roi_flank_size_bp))
            poi_by_chrom[chrom].append((end - roi_flank_size_bp, end + roi_flank_size_bp))

    for chrom in poi_by_chrom:
        intervals = sorted(poi_by_chrom[chrom])
        if not intervals: continue
        merged = []
        curr_start, curr_end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start <= curr_end:
                curr_end = max(curr_end, next_end)
            else:
                merged.append((curr_start, curr_end))
                curr_start, curr_end = next_start, next_end
        merged.append((curr_start, curr_end))
        poi_by_chrom[chrom] = merged
        print(f"DEBUG: bin_baf - Chromosome {chrom} has {len(merged)} ROI intervals", file=sys.stderr)

    def is_in_poi(chrom, pos):
        if chrom not in poi_by_chrom: return False
        for p_start, p_end in poi_by_chrom[chrom]:
            if pos >= p_start and pos <= p_end: return True
            if p_start > pos: break
        return False

    # Separate populations only for budget calculation
    roi_points = []
    normal_points = []
    for p in baf_list:
        if is_in_poi(p[0], p[1]): roi_points.append(p)
        else: normal_points.append(p)

    # Budget Decision
    roi_budget = target_data_points * roi_budget_fraction
    use_raw_roi = len(roi_points) <= roi_budget

    if use_raw_roi:
        target_normal = max(10, target_data_points - len(roi_points))
        k_normal = max(1, len(normal_points) / target_normal) if normal_points else 1
        k_roi = 1
    else:
        k_roi = (len(roi_points) + len(normal_points) / roi_resolution_factor) / target_data_points
        k_roi = max(1, k_roi)
        k_normal = k_roi * roi_resolution_factor

    print(f"DEBUG: [BAF] Total: {len(baf_list)}, ROI: {len(roi_points)}, Normal: {len(normal_points)}", file=sys.stderr)
    print(f"DEBUG: [BAF] Target: {target_data_points}, Budget: {roi_budget}, UseRaw: {use_raw_roi}", file=sys.stderr)
    print(f"DEBUG: [BAF] k_roi: {k_roi:.2f}, k_normal: {k_normal:.2f}", file=sys.stderr)

    binned_res = []
    current_bin = []
    current_bin_type = None

    def flush_bin(bin_to_flush):
        if not bin_to_flush: return
        n = len(bin_to_flush)
        start, end = bin_to_flush[0][1], bin_to_flush[-1][1]
        chrom = bin_to_flush[0][0]

        # Mode detection: Unified if noisy centered around 0.5, Split if bimodal
        # Histogram-like check for central peak [0.4, 0.6]
        mid_count = sum(1 for x in bin_to_flush if 0.4 <= x[2] <= 0.6)
        
        # If enough points are in the middle, keep it simple to avoid artificial splitting
        if mid_count > n * 0.4 or n < 4:
            mean = sum(x[2] for x in bin_to_flush) / n
            p = {"chromosome": chrom, "pos": (start + end) // 2, "start": start, "end": end, "baf": mean}
            if n > 1:
                p["mean"] = mean
                p["sd"] = math.sqrt(sum((x[2] - mean)**2 for x in bin_to_flush) / n)
            binned_res.append(p)
        else:
            # Bimodal (Split Hi/Lo)
            hi_sub = [x for x in bin_to_flush if x[2] >= 0.5]
            lo_sub = [x for x in bin_to_flush if x[2] < 0.5]
            for sub, tag in [(hi_sub, "Hi"), (lo_sub, "Lo")]:
                if not sub: continue
                sn = len(sub)
                s_mean = sum(x[2] for x in sub) / sn
                p = {"chromosome": chrom, "pos": (start + end) // 2, "start": start, "end": end, "baf": s_mean}
                if sn > 1:
                    p["mean"] = s_mean
                    p["sd"] = math.sqrt(sum((x[2] - s_mean)**2 for x in sub) / sn)
                binned_res.append(p)

    for chrom, pos, baf in baf_list:
        rtype = 'roi' if is_in_poi(chrom, pos) else 'normal'
        k_limit = k_roi if rtype == 'roi' else k_normal

        if current_bin:
            if rtype != current_bin_type or chrom != current_bin[0][0] or len(current_bin) >= k_limit:
                flush_bin(current_bin)
                current_bin = []
        
        if not current_bin:
            current_bin_type = rtype
        current_bin.append((chrom, pos, baf))

    flush_bin(current_bin)
    print(f"DEBUG: [BAF-Total] Result size: {len(binned_res)}", file=sys.stderr)
    return binned_res


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


def merge_cnv_dicts(
    dicts, baf, annotations, cytobands, chromosomes,
    filtered_cnvs, unfiltered_cnvs, gene_index=None, bin_params=None
):
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
            baf=[],
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

    if baf is not None:
        baf = list(baf)  # Consume generator
        # Combine all interesting regions for BAF binning
        # Annotations (genes) get full-region resolution
        # Segments (calls) get breakpoint resolution
        poi_regions = []
        for a in annotations:
            for item in a:
                poi_regions.append({"chromosome": normalize_chrom(item[0]), "start": item[1], "end": item[2], "full_region": True})
        for d in dicts:
            for s in d["segments"]:
                # Ensure chromosome is normalized from the loaded JSON
                s_normalized = {**s, "chromosome": normalize_chrom(s["chromosome"]), "full_region": False}
                poi_regions.append(s_normalized)

        # Bin ALL BAF data points globally using the budget
        global_binned_baf = bin_baf(baf, poi_regions, **(bin_params or {}))
        
        # Distribute binned points to chromosomes
        for b in global_binned_baf:
            chrom = b["chromosome"]
            if chrom in cnvs:
                # Remove chromosome key before adding to the per-chromosome list to keep JSON small
                point = {k:v for k,v in b.items() if k != "chromosome"}
                cnvs[chrom]["baf"].append(point)

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

    # Attach gene search index to first chromosome dictionary if available.
    # The frontend expects to find it in the first element of the chromosome array.
    if gene_index and cnvs:
        first_chrom = next(iter(cnvs))
        cnvs[first_chrom]["gene_search_index"] = gene_index

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
    baf = None
    if germline_vcf is not None:
        baf = get_baf(germline_vcf, skip_chromosomes)

    annot_parser = annotation_parser()
    annotations = []
    for filename in annotation_beds:
        annotations.append(list(annot_parser(filename, skip_chromosomes)))

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

    bin_params = {
        "roi_flank_size_bp": snakemake.params["roi_flank_size_bp"],
        "target_data_points": snakemake.params["target_data_points"],
        "roi_budget_fraction": snakemake.params["roi_budget_fraction"],
        "roi_resolution_factor": snakemake.params["roi_resolution_factor"],
    }

    cnvs = merge_cnv_dicts(
        cnv_dicts, baf, annotations, cytobands, fai, filtered_cnv_vcfs,
        unfiltered_cnv_vcfs, gene_index, bin_params
    )

    with open(output_file, "w") as f:
        print(json.dumps(cnvs, default=vars), file=f)


if __name__ == "__main__":
    main()
