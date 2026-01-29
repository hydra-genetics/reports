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

    # Sort ratios by position since we are not binning anymore
    ratios.sort(key=lambda x: (x["chromosome"], x["start"]))
    is_binned = False

    with open(output_filename, "w") as f:
        print(to_json(caller, ratios, segments, is_binned=is_binned), file=f)


if __name__ == "__main__":
    main()
