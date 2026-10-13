import re
import sys
import types
from pathlib import Path

import pandas as pd

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
TEMPLATE_DIR = TEST_DIR / "../../workflow/templates"
sys.path.insert(0, str(SCRIPT_DIR))
sys.modules.setdefault("snakemake", types.ModuleType("snakemake"))

import cnv_html_report  # noqa: E402
import general_html_report  # noqa: E402
import general_json_report  # noqa: E402

PURITY_PLOIDY = dict(tc_pathology=0.45, tc_purecn="0.39", tc_ichorcna=0.4635, ploidy_purecn=4.4389, ploidy_ichorcna=2.512)
UNITS = pd.DataFrame({"sample": ["S1"], "type": ["T"]})
PIPELINE = {"pipeline": {"version": "1.0.0"}}


def cnv_header(tmp_path, purity_ploidy, tc=0.39, tc_method="purecn"):
    json_file = tmp_path / "cnv.json"
    json_file.write_text("{}")

    class Wildcards:
        sample = "S1"
        type = "T"

    out = cnv_html_report.create_report(
        TEMPLATE_DIR / "cnv_html_report/index.html", json_file, [], [], True, [], tc, 4.44, tc_method, Wildcards(),
        False, False, False, False, False, purity_ploidy,
    )
    lines = re.findall(r'<p class="metadata-line">(.*?)</p>', out, re.S)
    return [re.sub(r"\s+", " ", line).strip() for line in lines]


def general_header(tc):
    config = dict(general_json_report.metadata("S1_T", PIPELINE, tc, UNITS, "ref.fa"), results=[])
    out = general_html_report.generate_report(
        str(TEMPLATE_DIR / "general_html_report/index.html"), config, 1, [], [], [], []
    )
    rows = re.findall(r"<td><b>(Tc [^<]+|Ploidy [^<]+)</b></td>\s*<td>([^<]*)</td>", out)
    return {name: value.strip() for name, value in rows}


def test_cnv_header_shows_all_estimates(tmp_path):
    lines = cnv_header(tmp_path, PURITY_PLOIDY)
    assert lines[0].endswith("Tumor cell content: 0.39 (purecn)")
    assert lines[1] == "TC pathology: 0.45 | TC PureCN: 0.39 | TC ichorCNA: 0.46 | Ploidy PureCN: 4.44 | Ploidy ichorCNA: 2.51"


def test_cnv_header_without_pathology_and_ichorcna(tmp_path):
    lines = cnv_header(tmp_path, dict(tc_pathology=None, tc_purecn="0.39", tc_ichorcna=None,
                                      ploidy_purecn=4.4389, ploidy_ichorcna=None))
    assert lines[1] == "TC PureCN: 0.39 | Ploidy PureCN: 4.44"


def test_cnv_header_keeps_zero_values(tmp_path):
    lines = cnv_header(tmp_path, dict(tc_pathology=0.0, tc_purecn=None, tc_ichorcna=0.0,
                                      ploidy_purecn=None, ploidy_ichorcna=2.15))
    assert lines[1] == "TC pathology: 0.00 | TC ichorCNA: 0.00 | Ploidy ichorCNA: 2.15"


def test_cnv_header_no_estimates(tmp_path):
    empty = dict(tc_pathology=None, tc_purecn=None, tc_ichorcna=None, ploidy_purecn=None, ploidy_ichorcna=None)
    assert len(cnv_header(tmp_path, empty)) == 1


def test_cnv_header_without_purity_ploidy(tmp_path):
    assert len(cnv_header(tmp_path, None)) == 1


def test_general_header_shows_all_estimates():
    assert general_header([0.45, "0.39", 4.4389, 0.4635, 2.512]) == {
        "Tc pathology": "0.45", "Tc purecn": "0.39", "Tc ichorCNA": "0.4635",
        "Ploidy purecn": "4.44", "Ploidy ichorCNA": "2.51",
    }


def test_general_header_hides_missing_new_values():
    assert general_header([0.45, "0.39", "NA", "NA", "NA"]) == {"Tc pathology": "0.45", "Tc purecn": "0.39"}


def test_general_header_accepts_two_values():
    assert general_header([0.45, "0.39"]) == {"Tc pathology": "0.45", "Tc purecn": "0.39"}


def test_general_header_keeps_zero_values():
    header = general_header([0.0, "NA", "NA", 0.0, 2.15])
    assert header["Tc pathology"] == "0.0"
    assert header["Tc ichorCNA"] == "0.0"
    assert header["Ploidy ichorCNA"] == "2.15"


def test_general_header_hidden_when_all_missing():
    assert general_header(["NA"] * 5) == {}


def test_cnv_header_without_tumor_content(tmp_path):
    lines = cnv_header(tmp_path, None, tc=None, tc_method="purecn")
    assert lines[0].endswith("Tumor cell content: N/A (purecn)")
