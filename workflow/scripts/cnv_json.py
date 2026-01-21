import collections
import csv
import functools
import json
from pathlib import Path
import sys


# The functions `parse_*_ratios` functions take a filename of a file containing
# copy number log2-ratios across the genome for a specific CNV caller. The
# functions `parse_*_segments` takes a filename of a file containing log2-ratio
# segments across the genome for a specific caller. The return value from both
# of these functions should be a list of dictionaries, where the dictionaries
# look like
#
# {
#     "chromosome": str,
#     "start": int,
#     "end": int,
#     "log2": float,
# }
#


PARSERS = collections.defaultdict(dict)


def normalize_chrom(chrom: str) -> str:
    """Ensure chromosome name starts with 'chr' exactly once."""
    c = str(chrom)
    if c.startswith("chr"):
        c = c[3:]
    return f"chr{c}"


def cnv_parser(file_format, header=True, skip=0, comment="#"):
    """
    Decorator for parsers of CNV result files. The first argument of
    the wrapped function should be a path to a file, and this argument
    is replaced with the contents of that file. How the content is represented
    depends on the file format:

    - tsv, csv: a generator over lines, each line being a list of values
    - vcf:

    arguments:
        file_format     a file path
        skip            the number of lines to skip before starting to read
                        the file
    """

    def decorator_cnv_parser(func):
        caller_filetype = func.__name__.split("_")[1:]
        caller = "_".join(caller_filetype[:-1])
        filetype = caller_filetype[-1]

        def line_generator(file, delim):
            for _ in range(skip):
                next(file)
            found_header = False
            for line in csv.reader(file, delimiter=delim):
                if line[0].strip()[0] == comment:
                    continue
                if header and not found_header:
                    found_header = True
                    continue
                yield line

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            f = None
            try:
                f = open(args[0], "r")
                if file_format == "tsv":
                    lines = line_generator(f, delim="\t")
                elif file_format == "csv":
                    lines = line_generator(f, delim=",")
                else:
                    raise IOError(f"invalid filetype: {file_format}")
                args = [lines] + list(args[1:])
                return func(*args, **kwargs)
            finally:
                if f is not None:
                    f.close()

        PARSERS[caller][filetype] = wrapper

    return decorator_cnv_parser


@cnv_parser("tsv", header=True)
def parse_cnvkit_ratios(file):
    ratios = []
    for line in file:
        ratios.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[5]),
            )
        )
    return ratios


@cnv_parser("tsv", header=True)
def parse_cnvkit_segments(file):
    segments = []
    for line in file:
        segments.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[4]),
            )
        )
    return segments


@cnv_parser("tsv", header=True, comment="@")
def parse_gatk_ratios(file):
    ratios = []
    for line in file:
        ratios.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[3]),
            )
        )
    return ratios


@cnv_parser("tsv", header=True, comment="@")
def parse_gatk_segments(file):
    segments = []
    for line in file:
        segments.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[4]),
            )
        )
    return segments


@cnv_parser("tsv", header=True)
def parse_jumble_ratios(file):
    ratios = []
    for line in file:
        ratios.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[5]),
            )
        )
    return ratios


@cnv_parser("tsv", header=True)
def parse_jumble_segments(file):
    segments = []
    for line in file:
        segments.append(
            dict(
                chromosome=normalize_chrom(line[0]),
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[4]),
            )
        )
    return segments


def bin_ratios(ratios, segments, roi_flank_size_bp=10000, target_data_points=50000, roi_budget_fraction=0.5, roi_resolution_factor=10):
    """
    Bin log2 ratios dynamically using a Budgeted ROI strategy.
    ROIs (genes) stay raw if they fit within a budget (X% of target).
    Otherwise, they are binned with higher resolution (Y times) than normal regions.
    """
    if not ratios:
        return []

    # Ensure ratios are sorted
    ratios.sort(key=lambda x: (x["chromosome"], x["start"]))

    # Group ROI by chromosome
    poi_by_chrom = collections.defaultdict(list)
    for s in segments:
        chrom = s["chromosome"]
        full_region = s.get("full_region", False)
        if full_region:
            poi_by_chrom[chrom].append((s["start"] - roi_flank_size_bp, s["end"] + roi_flank_size_bp))
        else:
            # Only breakpoints
            poi_by_chrom[chrom].append((s["start"] - roi_flank_size_bp, s["start"] + roi_flank_size_bp))
            poi_by_chrom[chrom].append((s["end"] - roi_flank_size_bp, s["end"] + roi_flank_size_bp))

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

    def is_in_poi(r):
        chrom = r["chromosome"]
        if chrom not in poi_by_chrom: return False
        r_start, r_end = r["start"], r["end"]
        for p_start, p_end in poi_by_chrom[chrom]:
            if r_start < p_end and r_end > p_start: return True
            if p_start > r_end: break
        return False

    # Separate populations
    roi_points = []
    normal_points = []
    for r in ratios:
        if is_in_poi(r):
            roi_points.append(r)
        else:
            normal_points.append(r)

    # Budget Decision
    roi_budget = target_data_points * roi_budget_fraction
    use_raw_roi = len(roi_points) <= roi_budget

    if use_raw_roi:
        target_normal = max(100, target_data_points - len(roi_points))
        k_normal = max(1, len(normal_points) / target_normal) if normal_points else 1
        k_roi = 1 # Keep raw
    else:
        k_roi = (len(roi_points) + len(normal_points) / roi_resolution_factor) / target_data_points
        k_roi = max(1, k_roi)
        k_normal = k_roi * roi_resolution_factor

    print(f"DEBUG: [log2] Total: {len(ratios)}, ROI: {len(roi_points)}, Normal: {len(normal_points)}", file=sys.stderr)
    print(f"DEBUG: [log2] Target: {target_data_points}, Budget: {roi_budget}, UseRaw: {use_raw_roi}", file=sys.stderr)
    print(f"DEBUG: [log2] k_roi: {k_roi:.2f}, k_normal: {k_normal:.2f}", file=sys.stderr)

    # Final Pass
    binned = []
    current_bin = []
    current_bin_type = None # 'roi' or 'normal'

    for r in ratios:
        in_roi = is_in_poi(r)
        rtype = 'roi' if in_roi else 'normal'
        k_limit = k_roi if in_roi else k_normal

        flush = False
        if current_bin:
            if rtype != current_bin_type or r["chromosome"] != current_bin[0]["chromosome"]:
                flush = True
            elif len(current_bin) >= k_limit:
                flush = True

        if flush and current_bin:
            n = len(current_bin)
            mean_log2 = sum(x["log2"] for x in current_bin) / n
            start, end = current_bin[0]["start"], current_bin[-1]["end"]
            binned.append({
                "chromosome": current_bin[0]["chromosome"],
                "start": start, "end": end, "pos": (start + end) // 2,
                "log2": mean_log2
            })
            current_bin = []

        if not current_bin:
            current_bin_type = rtype
        current_bin.append(r)

    if current_bin:
        n = len(current_bin)
        mean_log2 = sum(x["log2"] for x in current_bin) / n
        start, end = current_bin[0]["start"], current_bin[-1]["end"]
        binned.append({
            "chromosome": current_bin[0]["chromosome"],
            "start": start, "end": end, "pos": (start + end) // 2,
            "log2": mean_log2
        })

    return binned


def to_json(caller, ratios, segments, is_binned=False):
    json_dict = dict(
        caller=caller,
        ratios=ratios,
        segments=segments,
        is_binned=is_binned,
    )
    return json.dumps(json_dict)


def main():
    log = Path(snakemake.log[0])

    logfile = open(log, "w")
    sys.stdout = sys.stderr = logfile

    caller = snakemake.wildcards["caller"]
    ratio_filename = snakemake.input["ratios"]
    segment_filename = snakemake.input["segments"]
    annotation_filenames = snakemake.input["annotations"]

    output_filename = snakemake.output["json"]

    skip_chromosomes = snakemake.params["skip_chromosomes"]

    csv.field_size_limit(snakemake.params["csv_field_size_limit"])

    if caller not in PARSERS:
        print(f"error: no parser for {caller} implemented", file=sys.stderr)
        sys.exit(1)

    ratios = PARSERS[caller]["ratios"](ratio_filename)
    segments = PARSERS[caller]["segments"](segment_filename)

    # Parse annotations to include as ROIs
    annotations = []
    for annot_file in annotation_filenames:
        with open(annot_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    chrom = normalize_chrom(parts[0])
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        annotations.append({"chromosome": chrom, "start": start, "end": end})
                    except ValueError:
                        continue

    if skip_chromosomes is not None:
        ratios = [r for r in ratios if r["chromosome"] not in skip_chromosomes]
        segments = [s for s in segments if s["chromosome"] not in skip_chromosomes]
        annotations = [a for a in annotations if a["chromosome"] not in skip_chromosomes]

    # Perform smart binning for performance (especially for WGS)
    # Combine segments and annotations for POI identification
    # Annotations (genes) get full-region resolution
    # Segments (calls) get breakpoint resolution
    poi_regions = [
        {**s, "full_region": False} for s in segments
    ] + [
        {**a, "full_region": True} for a in annotations
    ]

    bin_params = {
        "roi_flank_size_bp": snakemake.params["roi_flank_size_bp"],
        "target_data_points": snakemake.params["target_data_points"],
        "roi_budget_fraction": snakemake.params["roi_budget_fraction"],
        "roi_resolution_factor": snakemake.params["roi_resolution_factor"],
    }

    original_ratio_count = len(ratios)
    ratios = bin_ratios(ratios, poi_regions, **bin_params)
    is_binned = len(ratios) < original_ratio_count

    with open(output_filename, "w") as f:
        print(to_json(caller, ratios, segments, is_binned=is_binned), file=f)


if __name__ == "__main__":
    main()
