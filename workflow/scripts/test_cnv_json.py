from pathlib import Path
import pytest
import sys

SCRIPT_DIR = Path(__file__).parent.resolve()
INTEGRATION_DIR = SCRIPT_DIR / "../../.tests/integration"
sys.path.insert(0, str(SCRIPT_DIR))

import cnv_json  # noqa


@pytest.fixture
def cnvkit_segment_file():
    return INTEGRATION_DIR / "cnv_sv/cnvkit_batch/sample1/sample1_T.cns"


def test_existing_parsers():
    assert "cnvkit" in cnv_json.PARSERS
    assert "gatk" in cnv_json.PARSERS


def test_parse_cnvkit_segments(cnvkit_segment_file):
    segments = cnv_json.PARSERS["cnvkit"]["segments"](cnvkit_segment_file)
    assert len(segments) == 13
    assert all(x in segments[0] for x in ["chromosome", "start", "end", "log2"])
