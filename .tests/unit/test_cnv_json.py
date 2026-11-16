from pathlib import Path
import pytest
import sys

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
INTEGRATION_DIR = TEST_DIR / "../integration"
sys.path.insert(0, str(SCRIPT_DIR))

import cnv_json  # noqa


@pytest.fixture
def cnvkit_segment_file():
    return INTEGRATION_DIR / "cnv_sv/cnvkit_batch/sample1/sample1_T.cns"


@pytest.fixture
def jumble_ratio_file():
    return INTEGRATION_DIR / "cnv_sv/jumble_run/sample1_T/sample1_T.cnr"


@pytest.fixture
def jumble_segment_file():
    return INTEGRATION_DIR / "cnv_sv/jumble_run/sample1_T/sample1_T.cns"


@pytest.fixture
def gatk_ratio_file():
    return INTEGRATION_DIR / "cnv_sv/gatk_denoise_read_counts/sample1_T.clean.denoisedCR.tsv"


@pytest.fixture
def gatk_segment_file():
    return INTEGRATION_DIR / "cnv_sv/gatk_model_segments/sample1_T.clean.cr.seg"


def test_existing_parsers():
    assert "cnvkit" in cnv_json.PARSERS
    assert "gatk" in cnv_json.PARSERS
    assert "jumble" in cnv_json.PARSERS


def test_parse_cnvkit_segments(cnvkit_segment_file):
    segments = cnv_json.PARSERS["cnvkit"]["segments"](cnvkit_segment_file)
    assert len(segments) == 13
    assert all(x in segments[0] for x in ["chromosome", "start", "end", "log2"])


def test_parse_jumble_ratios(jumble_ratio_file):
    # sample1_T.cnr has the current Jumble .cnr header (with bed_name/gc/count/type,
    # which don't exist in cnvkit's/older Jumble's .cnr) - column-name-based parsing
    # must still find log2 correctly regardless of its position in the header.
    ratios = cnv_json.PARSERS["jumble"]["ratios"](jumble_ratio_file)
    assert len(ratios) == 11727
    assert all(x in ratios[0] for x in ["chromosome", "start", "end", "log2"])
    assert ratios[0]["chromosome"] == "chrA"
    assert ratios[0]["log2"] == pytest.approx(1.1927529757296118)


def test_parse_jumble_segments(jumble_segment_file):
    # sample1_T.cns has the current Jumble .cns header (with band/relevance instead
    # of cnvkit's weight/ci_lo/ci_hi) - same column-name-based parsing requirement.
    segments = cnv_json.PARSERS["jumble"]["segments"](jumble_segment_file)
    assert len(segments) == 13
    assert all(x in segments[0] for x in ["chromosome", "start", "end", "log2"])
    assert segments[1]["log2"] == pytest.approx(-3.530150424344391)


def test_parse_gatk_ratios(gatk_ratio_file):
    # GATK's own header uses CONTIG/LOG2_COPY_RATIO rather than
    # chromosome/log2 - the parser must map those, not assume the same names
    # cnvkit/jumble use.
    ratios = cnv_json.PARSERS["gatk"]["ratios"](gatk_ratio_file)
    assert len(ratios) == 600
    assert all(x in ratios[0] for x in ["chromosome", "start", "end", "log2"])
    assert ratios[0]["chromosome"] == "chrA"
    assert ratios[0]["log2"] == pytest.approx(0.28)


def test_parse_gatk_segments(gatk_segment_file):
    segments = cnv_json.PARSERS["gatk"]["segments"](gatk_segment_file)
    assert len(segments) == 13
    assert all(x in segments[0] for x in ["chromosome", "start", "end", "log2"])
    assert segments[0]["log2"] == pytest.approx(0.77)
