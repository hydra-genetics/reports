import pytest
from pathlib import Path
import sys

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from cnv_html_report import get_sample_name


def test_get_sample_name():
    class MockWildcards:
        def __init__(self, sample):
            self.sample = sample

    assert get_sample_name(MockWildcards("sample")) == "sample"
    assert get_sample_name(MockWildcards("sample.type")) == "sample.type"
    assert get_sample_name(MockWildcards("2021.CC.13_T")) == "2021.CC.13_T"
