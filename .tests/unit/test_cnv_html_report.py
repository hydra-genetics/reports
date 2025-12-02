import pytest
from pathlib import Path
import sys

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from cnv_html_report import get_sample_name


def test_get_sample_name():
    # simple case
    assert get_sample_name("sample.json") == "sample"
    # nested suffixes
    assert get_sample_name("sample.type.json") == "sample.type"
    assert get_sample_name("sample.type.json.gz") == "sample.type"
    # complex sample name
    assert get_sample_name("reports/cnv_html_report/2021.CC.13_T.tc_method.merged.json") == "2021.CC.13_T"
