from dataclasses import dataclass
import os
import sys
from typing import Dict, List
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_cnv_json import CNV, filter_chr_cnvs  # noqa


class TestMergeCnvJson(unittest.TestCase):
    def test_cnv_filter(self):
        """
        The reason for this test was an issue where the filter for
        the results table behaved weirdly in that segments that
        should have shown up in the filtered table did not.
        """

        @dataclass
        class TestCase:
            name: str
            input: Dict[str, List[CNV]]
            expected_names: List[str]
            expected_lens: Dict[str, int]
            expected_filter: Dict[str, List[bool]]

        testcases = [
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
            TestCase(
                name="table filter bug",
                input={
                    "cnvkit": [
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene1",
                            start=1000,
                            length=2500,
                            type="COPY_NORMAL",
                            copy_number=2,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene2",
                            start=4000,
                            length=1000,
                            type="DUP",
                            copy_number=9,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene3",
                            start=5800,
                            length=400,
                            type="COPY_NORMAL",
                            copy_number=2,
                            baf=0.6,
                        ),
                    ],
                    "gatk": [
                        CNV(
                            caller="gatk",
                            chromosome="chr1",
                            genes=["gene1", "gene2", "gene3"],
                            start=10,
                            length=7190,
                            type="COPY_NORMAL",
                            copy_number=2,
                            baf=0.6,
                        )
                    ],
                },
                expected_lens={"cnvkit": 3, "gatk": 1},
                expected_names=["cnvkit", "gatk"],
                expected_filter={"cnvkit": [False, True, False], "gatk": [True]},
            ),
            TestCase(
                name="example from real report",
                input={
                    "cnvkit": [
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene1"],
                            start=132_309,
                            length=3_598_391,
                            type="DEL",
                            copy_number=1.41,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene2"],
                            start=3_731_199,
                            length=31_683_311,
                            type="COPY_NORMAL",
                            copy_number=1.85,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene3"],
                            start=66_483_605,
                            length=4_031_951,
                            type="DEL",
                            copy_number=1.28,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene4"],
                            start=87_570_960,
                            length=2_524_518,
                            type="DEL",
                            copy_number=1.20,
                            baf=0.5,
                        ),
                    ],
                    "gatk": [
                        CNV(
                            caller="gatk",
                            chromosome="chr16",
                            genes=["gene1", "gene2", "gene3", "gene4"],
                            start=132_060,
                            length=89_963_668,
                            type="COPY_NORMAL",
                            copy_number=0.99,
                            baf=0.5,
                        )
                    ],
                },
                expected_names=["cnvkit", "gatk"],
                expected_lens={"cnvkit": 4, "gatk": 1},
                expected_filter={"cnvkit": [True, False, True, True], "gatk": [True]},
            ),
        ]

        for case in testcases:
            filtered_cnvs = {}
            for caller, cnvs in case.input.items():
                f = list(filter(lambda x: x.type != "COPY_NORMAL", cnvs))
                if len(f) > 0:
                    filtered_cnvs[caller] = f

            m = filter_chr_cnvs(
                unfiltered_cnvs=case.input,
                filtered_cnvs=filtered_cnvs,
            )

            assert list(m.keys()) == case.expected_names

            for caller, cnvs in m.items():
                assert len(cnvs) == case.expected_lens[caller]
                for cnv, expected_filter in zip(cnvs, case.expected_filter[caller]):
                    assert cnv["passed_filter"] == expected_filter
