import pytest # type: ignore
from pathlib import Path
import sys
import pandas as pd


TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
RESULTS_DIR = TEST_DIR / "../../results"
sys.path.insert(0, str(SCRIPT_DIR))

import pysam

from compile_xlsx_report import ( # type: ignore
    parse_info,
    parse_vcf_line,
    parse_format,
    parse_cnvkit_vcf_line,
    parse_version_from_container,
    get_genes_from_bed,
    process_sv_vcf,
)


# --- Tests for parse_info ---
@pytest.mark.parametrize(
    "info_str, expected",
    [
        # Normal case
        ("FAU=46;FCU=28", {"FAU": "46", "FCU": "28"}),
        # empty input
        ("", {}),
        # malformed input
        ("FAU;FCU=28", {"FCU": "28"}),
    ],
)
def test_parse_info(info_str, expected):
    assert parse_info(info_str) == expected


# --- Tests for parse_format ---
@pytest.mark.parametrize(
    "format_str, sample_str, expected",
    [
        # Normal case
        ("GT:GQ", "0/1:76", {"GT": "0/1", "GQ": "76"}),
        # Short input strings (one item in each)
        ("GT", "0/1", {"GT": "0/1"}),
        # Both empty strings
        ("", "", {"": ""}),
        # Non-string inputs → should raise AttributeError because .split() fails
        (None, "0/1:76", pytest.raises(AttributeError)),
        ("GT:GQ", None, pytest.raises(AttributeError)),
        (123, "0/1:76", pytest.raises(AttributeError)),
        ("GT:GQ", 456, pytest.raises(AttributeError)),
        (None, None, pytest.raises(AttributeError)),
        (123, 123, pytest.raises(AttributeError)),
        # Unequal length of input strings
        ("GT:GQ:BD", "0/1:76", pytest.raises(ValueError)),
        # one input is empty the other is not
        ("", "0/1:76", pytest.raises(ValueError)),
        ("GT:GQ:BD", "", pytest.raises(ValueError)),
    ],
)
def test_parse_format(format_str, sample_str, expected):
    if isinstance(expected, dict):
        assert parse_format(format_str, sample_str) == expected
    else:
        with expected:
            parse_format(format_str, sample_str)


# --- Tests for parse_vcf_line ---
@pytest.mark.parametrize(
    "line, vep_fields, format_fields, expected",
    [
        # Normal case — no CALLER in INFO, CALLER defaults to empty string
        (
            "chr1\t12345\t.\tA\tT\t33\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP"],
            {
                "CHROM": "chr1",
                "POS": "12345",
                "REF": "A",
                "ALT": "T",
                "QUAL": "33",
                "FILTER": "PASS",
                "CALLER": "",
                "GT": "0/1",
                "DP": "35",
                "Gene": "gene1",
                "Impact": "impact1",
            },
        ),
        # CALLER tag present in INFO
        (
            "chr1\t12345\t.\tA\tT\t33\tPASS\tCALLER=clairs_to;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP"],
            {
                "CHROM": "chr1",
                "POS": "12345",
                "REF": "A",
                "ALT": "T",
                "QUAL": "33",
                "FILTER": "PASS",
                "CALLER": "clairs_to",
                "GT": "0/1",
                "DP": "35",
                "Gene": "gene1",
                "Impact": "impact1",
            },
        ),
        # CALLER=deepsomatic in INFO
        (
            "chr1\t200\t.\tG\tC\t50\tPASS\tCALLER=deepsomatic;CSQ=geneX|HIGH\tGT:DP\t1/1:40",
            ["Gene", "Impact"],
            ["GT", "DP"],
            {
                "CHROM": "chr1",
                "POS": "200",
                "REF": "G",
                "ALT": "C",
                "QUAL": "50",
                "FILTER": "PASS",
                "CALLER": "deepsomatic",
                "GT": "1/1",
                "DP": "40",
                "Gene": "geneX",
                "Impact": "HIGH",
            },
        ),
        # Too few columns in VCF line
        ("chr1\t12345\t.\tA\tT", [], [], pytest.raises(ValueError)),
        # GQ missing in the FORMAT string — parsed as empty string rather than error
        # (e.g. DeepSomatic records lack AF; missing fields must not abort parsing)
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:25",
            ["Gene", "Impact"],
            ["GT", "DP", "GQ"],
            {
                "CHROM": "chr1",
                "POS": "12345",
                "REF": "A",
                "ALT": "T",
                "QUAL": ".",
                "FILTER": "PASS",
                "CALLER": "",
                "Gene": "gene1",
                "Impact": "impact1",
                "GT": "0/1",
                "DP": "25",
                "GQ": "",
            },
        ),
        # You request too many VEP fields from the VCF file
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact", "Substitution"],  # Substitution is missing
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
        # CSQ field is missing in VCF line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;\tGT:DP\t0/1:45",
            ["Gene", "Impact"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_vcf_line(line, vep_fields, format_fields, expected):
    if isinstance(expected, dict):
        assert (
            parse_vcf_line(line, vep_fields, format_fields)
            == expected
        )
    else:
        with expected:
            parse_vcf_line(line, vep_fields, format_fields)



@pytest.mark.parametrize(
    "line, expected",
    [
        # Normal case
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            {
                "CHROM": "chr1",
                "POS": 10,
                "VARIANT_TYPE": "COPY_NORMAL",
                "GENE": "gene1,gene2",
                "END": 125,
                "SVLEN": 124,
                "LOG_ODDS_RATIO": -0.1,
                "CORR_CN": 2.0,
                "PROBES": 2,
                "BAF": 0.3,
                "GT": "0/0",
                "CNQ": 80,
                "DP": 3.5,
            },
        ),
        # Genes in INFO is missing (it's OK when the variant is DEL)
        (
            "chr1\t10\t.\tN\t<DEL>\t.\t.\tSVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            {
                "CHROM": "chr1",
                "POS": 10,
                "VARIANT_TYPE": "DEL",
                # Genes missing → empty string
                "GENE": "",
                "END": 125,
                "SVLEN": 124,
                "LOG_ODDS_RATIO": -0.1,
                "CORR_CN": 2.0,
                "PROBES": 2,
                "BAF": 0.3,
                "GT": "0/0",
                "CNQ": 80,
                "DP": 3.5,
            },
        ),
        # Line is too short
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # NO info field
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\t?\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # One of the expected values in INFO is missing (SVLEN, SVTYPE etc)
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # One of the expected values in FORMAT is missing (GT, GQ etc)
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CNQ\t0/0:2.0",
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_cnvkit_vcf_line(line, expected):
    if isinstance(expected, dict):
        assert parse_cnvkit_vcf_line(line) == expected
    else:
        with expected:
            parse_cnvkit_vcf_line(line)


# --- Tests for parse_version_from_container ---
@pytest.mark.parametrize(
    "container, expected",
    [
        ("docker://hydragenetics/severus:1.5", "1.5"),
        ("docker://hydragenetics/pbmm2:1.16", "1.16"),
        ("docker://hkubal/clairs-to:latest", "latest"),
        ("docker://org/tool:", ""),
        ("someimage", "unknown"),
        ("", "unknown"),
    ],
)
def test_parse_version_from_container(container, expected):
    assert parse_version_from_container(container) == expected


# --- Integration test for process_sv_vcf ---
COLO829_VCF = RESULTS_DIR / "COLO829_T.vcf"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_columns():
    """process_sv_vcf must return all expected columns including the four population frequency columns."""
    df = process_sv_vcf(str(COLO829_VCF))
    expected_cols = [
        "CHROM", "POS", "END", "TYPE", "SVLEN", "ALT", "FILTER", "CALLER",
        "COVERAGE", "SUPPORT", "STRAND", "VAF", "GENOTYPE", "GENOME QUALITY",
        "DEPTH REF", "DEPTH TRANS",
        "GNOMAD_AC", "GNOMAD_AF", "CUSTOM_AC", "CUSTOM_AF",
    ]
    assert list(df.columns) == expected_cols
    assert len(df) > 0


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_types():
    """Population frequency columns must be numeric (not string) where present.
    pandas coerces int+None columns to float64, which is fine for Excel rendering."""
    df = process_sv_vcf(str(COLO829_VCF))
    for col in ["GNOMAD_AC", "GNOMAD_AF", "CUSTOM_AC", "CUSTOM_AF"]:
        assert pd.api.types.is_numeric_dtype(df[col]), \
            f"Column {col} should be numeric, got dtype {df[col].dtype}"
    # Spot-check: at least one row has a non-null GNOMAD_AC value
    assert df["GNOMAD_AC"].notna().any(), "Expected at least one non-null GNOMAD_AC"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_sample_selection():
    """
    Severus records should have depth from the haplotagged sample (DR/DV present).
    PBSV records should have depth from the AD fallback.
    """
    df = process_sv_vcf(str(COLO829_VCF))
    severus_rows = df[df["CALLER"] == "severus"]
    pbsv_rows = df[df["CALLER"] == "pbsv"]
    assert len(severus_rows) > 0, "Expected at least one severus record in fixture"
    assert len(pbsv_rows) > 0, "Expected at least one pbsv record in fixture"
    # Severus: DR and DV are explicit FORMAT fields → non-empty
    assert (severus_rows["DEPTH TRANS"] != "").all(), "Severus DEPTH TRANS should be non-empty"
    # PBSV: depth comes from AD fallback → non-empty
    assert (pbsv_rows["DEPTH TRANS"] != "").all(), "PBSV DEPTH TRANS should be non-empty (AD fallback)"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_severus_strand_vaf():
    """Severus records should have STRANDS → STRAND populated and VAF from FORMAT."""
    df = process_sv_vcf(str(COLO829_VCF))
    severus_rows = df[df["CALLER"] == "severus"]
    assert len(severus_rows) > 0, "Expected at least one severus record in fixture"
    assert (severus_rows["STRAND"] != "").all(), "Severus STRAND should be non-empty (from STRANDS field)"
    assert (severus_rows["VAF"] != "").all(), "Severus VAF should be non-empty (from FORMAT)"


# --- Tests for get_genes_from_bed ---
def test_get_genes_from_bed(tmp_path):
    # Test valid BED file
    bed_content = (
        "chr1\t100\t200\tGENE1\n"
        "chr2\t300\t400\tGENE2\tignore\n"
        "#comment line\n"
        "\n"
        "chr3\t500\t600\tGENE3\n"
    )
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(bed_content)
    
    genes = get_genes_from_bed(str(bed_file))
    assert genes == {"GENE1", "GENE2", "GENE3"}

    # Test missing file
    assert get_genes_from_bed("non_existent.bed") == set()

    # Test empty path
    assert get_genes_from_bed("") == set()

    # Test malformed line (fewer than 4 columns)
    malformed_content = "chr1\t100\t200\n"
    malformed_file = tmp_path / "malformed.bed"
    malformed_file.write_text(malformed_content)
    assert get_genes_from_bed(str(malformed_file)) == set()


def test_parse_vcf_line_caller_missing_defaults_to_empty():
    """When CALLER is absent from INFO, parse_vcf_line returns CALLER=''."""
    line = "chr1\t500\t.\tC\tG\t.\tPASS\tCSQ=geneA|LOW\tGT:DP\t0/1:20"
    result = parse_vcf_line(line, ["Gene", "Impact"], ["GT", "DP"])
    assert "CALLER" in result
    assert result["CALLER"] == ""


def test_parse_vcf_line_caller_present():
    """When CALLER=deepsomatic is in INFO, parse_vcf_line captures it correctly."""
    line = "chr7\t100\t.\tA\tT\t.\tPASS\tCALLER=deepsomatic;CSQ=tp53|HIGH\tGT:DP\t0/1:30"
    result = parse_vcf_line(line, ["Gene", "Impact"], ["GT", "DP"])
    assert result["CALLER"] == "deepsomatic"


# --- Regression test: SVLEN tuple unwrapping ---

_SV_VCF_TUPLE_SVLEN = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr14,length=107043718>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length (dot=tuple in pysam)">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strands">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=svdb_origin,Number=1,Type=String,Description="Caller origin">
##INFO=<ID=SUPPORT,Number=.,Type=Integer,Description="Support">
##INFO=<ID=VAF,Number=.,Type=Float,Description="VAF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Ref depth">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Alt depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr14\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-5000;END=6000;svdb_origin=pbsv;SUPPORT=10;VAF=0.5\tGT:DR:DV\t0/1:20:10
"""


def test_process_sv_vcf_svlen_tuple_parsed_as_number(tmp_path):
    """SVLEN declared as Number=. causes pysam to return a tuple.
    process_sv_vcf must unwrap it so pd.to_numeric succeeds and the record
    is not silently dropped by the SVLEN >= sv_min_len filter."""
    vcf_path = tmp_path / "sv_tuple_svlen.vcf"
    vcf_path.write_text(_SV_VCF_TUPLE_SVLEN)
    gz_path = str(tmp_path / "sv_tuple_svlen.vcf.gz")
    pysam.tabix_compress(str(vcf_path), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)

    df = process_sv_vcf(gz_path)
    assert len(df) == 1, "Expected one record in the output DataFrame"

    svlen_numeric = pd.to_numeric(df["SVLEN"], errors="coerce")
    assert svlen_numeric.notna().all(), (
        f"SVLEN could not be parsed as a number — raw value was {df['SVLEN'].iloc[0]!r}; "
        "likely a pysam tuple was not unwrapped before str() conversion"
    )
    assert abs(svlen_numeric.iloc[0]) == 5000

    # SUPPORT and VAF should also be unwrapped
    assert df["SUPPORT"].iloc[0] == "10", f"SUPPORT was {df['SUPPORT'].iloc[0]!r}, expected '10'"
    assert df["VAF"].iloc[0] == "0.50", f"VAF was {df['VAF'].iloc[0]!r}, expected '0.50'"
    assert df["END"].iloc[0] == "6000", f"END was {df['END'].iloc[0]!r}, expected '6000'"


# --- Regression test: PBSV VAF derived from AD ---

_SV_VCF_PBSV_NO_VAF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr14,length=107043718>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=svdb_origin,Number=1,Type=String,Description="Caller origin">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr14\t5000\tpbsv.BND.1\tN\t<BND>\t.\tPASS\tSVTYPE=BND;SVLEN=10000;END=15000;svdb_origin=pbsv\tGT:AD:DP\t0/1:30,10:40
"""


def test_process_sv_vcf_pbsv_vaf_derived_from_ad(tmp_path):
    """PBSV records carry FORMAT/AD (ref, alt) but no FORMAT/VAF.
    process_sv_vcf must derive VAF from AD so the SV and Translocations tabs
    are not silently left with empty VAF when only Severus+PBSV are merged."""
    vcf_path = tmp_path / "sv_pbsv_no_vaf.vcf"
    vcf_path.write_text(_SV_VCF_PBSV_NO_VAF)
    gz_path = str(tmp_path / "sv_pbsv_no_vaf.vcf.gz")
    pysam.tabix_compress(str(vcf_path), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)

    df = process_sv_vcf(gz_path)
    assert len(df) == 1, "Expected one record"

    vaf_str = df["VAF"].iloc[0]
    assert vaf_str != "", "VAF must not be empty for a PBSV record with AD"
    assert vaf_str == "0.25", (
        f"Expected VAF = alt_depth/total_depth formatted to 2 d.p. = '0.25', got {vaf_str!r}"
    )
    assert df["SUPPORT"].iloc[0] == "10", (
        f"Expected SUPPORT = AD[1] = 10 for a PBSV record, got {df['SUPPORT'].iloc[0]!r}"
    )
