from dataclasses import dataclass
import os
import sys
from typing import Any, Dict, List, Tuple
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_cnv_json import CNV, filter_chr_cnvs, merge_cnv_dicts  # noqa


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

    def test_merge_cnv_dicts(self):
        @dataclass
        class TestCase:
            chromosomes: List[Tuple[str, int]]
            callers: List[Dict[str, str]]
            unfiltered_cnvs: List[List[CNV]]
            filtered_cnvs: List[List[CNV]]
            expected_cnvs: Dict[str, Any]

        testcases = [
            TestCase(
                chromosomes=[("chr1", 1000)],
                callers=[
                    {
                        "caller": "cnvkit",
                        "segments": [],
                        "ratios": [],
                    },
                ],
                unfiltered_cnvs=[
                    # File with amplifications
                    {
                        "chr1": {
                            "cnvkit": [
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene1"],
                                    start=100,
                                    length=100,
                                    type="DUP",
                                    copy_number=6,
                                    baf=0.3,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    copy_number=2,
                                    baf=0.5,
                                ),
                            ],
                        },
                    },
                    # File with deletions
                    {
                        "chr1": {
                            "cnvkit": [
                                # This variant is present in the other file too,
                                # but we don't want duplicates in the final data.
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    copy_number=2,
                                    baf=0.5,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene1"],
                                    start=200,
                                    length=100,
                                    type="COPY_NORMAL",
                                    copy_number=1,
                                    baf=0.5,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=600,
                                    length=100,
                                    type="DEL",
                                    copy_number=1,
                                    baf=0.3,
                                ),
                            ],
                        },
                    },
                ],
                filtered_cnvs=[
                    # Filtered amplifications
                    {
                        "chr1": {
                            "cnvkit": [
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene1"],
                                    start=100,
                                    length=100,
                                    type="DUP",
                                    copy_number=6,
                                    baf=0.3,
                                ),
                            ],
                        },
                    },
                    # Filtered deletions
                    {
                        "chr1": {
                            "cnvkit": [
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=600,
                                    length=100,
                                    type="DEL",
                                    copy_number=1,
                                    baf=0.3,
                                ),
                            ],
                        },
                    },
                ],
                expected_cnvs={
                    "chr1": {
                        "cnvkit": {
                            # All four calls should end up in the results...
                            "all_count": 4,
                            # ... and two of them should pass the filtering
                            "pass_filter_count": 2,
                        }
                    }
                },
            )
        ]

        for case in testcases:
            d = merge_cnv_dicts(case.callers, [], [], [], case.chromosomes, case.filtered_cnvs, case.unfiltered_cnvs)

            for chromosome in case.expected_cnvs:
                chromosomedata = [x for x in d if x["chromosome"] == chromosome]
                assert len(chromosomedata) == 1

                for caller in case.expected_cnvs[chromosome]:
                    callerdata = [x for x in chromosomedata[0]["callers"] if x["name"] == caller]
                    assert len(callerdata) == 1

                    assert case.expected_cnvs[chromosome][caller]["all_count"] == len(callerdata[0]["cnvs"])

                    n_pass_filter = sum(x["passed_filter"] for x in callerdata[0]["cnvs"])
                    assert n_pass_filter == case.expected_cnvs[chromosome][caller]["pass_filter_count"]
