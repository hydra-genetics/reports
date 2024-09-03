import json
import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_cnv_json import CNV, merge_cnv_dicts  # noqa


class TestTableFilter(unittest.TestCase):
    def test_multiple_overlapping_segments(self):
        lr_dicts = [
            {
                "caller": "cnvkit",
                "ratios": [],
                "segments": [],
            },
            {
                "caller": "gatk",
                "ratios": [],
                "segments": [],
            },
        ]
        # These are the segments that the test is simulating.
        #
        # In this case the first cnvkit segment should not pass the
        # filter. The second cnvkit segment should pass since it's a
        # duplication, and the GATK segment should also pass since it
        # overlaps the cnvkit segment that does pass.
        #
        # n: copy normal, ^: duplication
        #
        # cnvkit:           nnnnnnnnnn          ^^^^^^^^^^
        # gatk:   nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
        # chr1:   __________________________________________________
        unfiltered_cnvs = [
            {
                "chr1": {
                    "cnvkit": [
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene1",
                            start=1000,
                            length=1000,
                            type="COPY_NORMAL",
                            copy_number=2,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene3",
                            start=3000,
                            length=1000,
                            type="DUP",
                            copy_number=9,
                            baf=0.95,
                        ),
                    ],
                    "gatk": [
                        CNV(
                            caller="gatk",
                            chromosome="chr1",
                            genes=["gene1", "gene2", "gene3"],
                            start=10,
                            length=3500,
                            type="COPY_NORMAL",
                            copy_number=2,
                            baf=0.5,
                        )
                    ],
                }
            }
        ]
        filtered_cnvs = [
            {
                "chr1": {
                    "cnvkit": list(filter(lambda x: x.type != "COPY_NORMAL", unfiltered_cnvs[0]["chr1"]["cnvkit"])),
                    "gatk": list(filter(lambda x: x.type != "COPY_NORMAL", unfiltered_cnvs[0]["chr1"]["gatk"])),
                }
            }
        ]
        chromosomes = [("chr1", 5000)]

        m = merge_cnv_dicts(
            dicts=lr_dicts,
            vaf=None,
            annotations=[],
            cytobands=[],
            chromosomes=chromosomes,
            filtered_cnvs=filtered_cnvs,
            unfiltered_cnvs=unfiltered_cnvs,
        )

        # Ordered dicts are guaranteed for python 3.7+
        assert m[0]["callers"][0]["name"] == "cnvkit"
        assert m[0]["callers"][1]["name"] == "gatk"

        # cnvkit segments
        assert len(m[0]["callers"][0]["cnvs"]) == 2, "expected 2 cnvkit segments"
        assert m[0]["callers"][0]["cnvs"][0]["passed_filter"] is False, "first cnvkit segment should fail filter"
        assert m[0]["callers"][0]["cnvs"][1]["passed_filter"] is True, "second cnvkit segment should pass filter"

        # gatk segment
        assert len(m[0]["callers"][1]["cnvs"]) == 1, "expected 1 gatk segment"
        assert m[0]["callers"][1]["cnvs"][0]["passed_filter"] is True, "gatk segment should pass filter"
