import pytest
from workflow.scripts.cnv_html_report import get_sample_name


def test_get_sample_name():
    # simple case
    assert get_sample_name("test.vcf") == "test"
    # nested suffixes
    assert get_sample_name("test.vcf.gz") == "test"
    assert get_sample_name("test.vcf.gz.gz") == "test"
    # complex sample name
    assert get_sample_name("reports/cnv_html_report/2021.CC.13_T.tc_method.merged.json") == "2021.CC.13_T"