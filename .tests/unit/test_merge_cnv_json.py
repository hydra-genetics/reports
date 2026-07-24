from collections import defaultdict
from dataclasses import dataclass
import os
import sys
import tempfile
from typing import Any, Dict, List, Tuple
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_cnv_json import (  # noqa
    CNV,
    filter_chr_cnvs,
    merge_cnv_dicts,
    get_cnvs,
    merge_cnv_calls,
    sort_cnvs,
    evaluate_filter_condition,
    passes_table_filter,
    classify_cnv,
)


class TestCNV(unittest.TestCase):
    def test_membership(self):
        cnv = CNV(
            caller="cnvkit",
            chromosome="chr1",
            start=100,
            length=100,
            type="COPY_NORMAL",
            genes=[],
            cn=2,
            baf=0.5,
            passed_filter=True,
        )

        cnv_collection = [
            CNV(
                caller="cnvkit",
                chromosome="chr1",
                start=100,
                length=100,
                type="COPY_NORMAL",
                genes=[],
                cn=2,
                baf=0.5,
            ),
        ]

        assert cnv in cnv_collection


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
            # cnvkit:          nnnnnnnnnnnnnnnnnnnnnnnnn     ^^^^^^^^^^^       nnnn
            # gatk:    nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
            # chr1:   ________________________________________________________________________
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
                            cn=2,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene2",
                            start=4000,
                            length=1000,
                            type="DUP",
                            cn=9,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr1",
                            genes="gene3",
                            start=5800,
                            length=400,
                            type="COPY_NORMAL",
                            cn=2,
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
                            cn=2,
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
                            cn=1.41,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene2"],
                            start=3_731_199,
                            length=31_683_311,
                            type="COPY_NORMAL",
                            cn=1.85,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene3"],
                            start=66_483_605,
                            length=4_031_951,
                            type="DEL",
                            cn=1.28,
                            baf=0.5,
                        ),
                        CNV(
                            caller="cnvkit",
                            chromosome="chr16",
                            genes=["gene4"],
                            start=87_570_960,
                            length=2_524_518,
                            type="DEL",
                            cn=1.20,
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
                            cn=0.99,
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
                    assert cnv.passed_filter == expected_filter

    def test_get_cnvs_baf_zero_is_missing_not_a_real_measurement(self):
        # Regression test for two real-world cases: a genuinely measured,
        # strongly-skewed BAF can be negative (e.g. from a misestimated TC)
        # and must still clamp to a real 0.0 — but a raw BAF of EXACTLY 0
        # (Jumble's permanent placeholder, or a CNVkit call with too few
        # PROBES to compute a meaningful BAF) must become None/missing
        # instead, so it isn't mistaken for a real "fully skewed" measurement.
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr3>
##contig=<ID=chr8>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="x">
##INFO=<ID=END,Number=1,Type=Integer,Description="x">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="x">
##INFO=<ID=CALLER,Number=1,Type=String,Description="x">
##INFO=<ID=CORR_CN,Number=1,Type=Float,Description="x">
##INFO=<ID=BAF,Number=1,Type=Float,Description="x">
##INFO=<ID=Genes,Number=.,Type=String,Description="x">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr3\t10077950\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=GENEX;SVTYPE=COPY_NORMAL;END=10191669;SVLEN=113720;CALLER=cnvkit;CORR_CN=2.21;BAF=-0.448501
chr8\t43475717\t.\tN\t<DUP>\t.\t.\tGenes=GENEY;SVTYPE=DUP;END=43930141;SVLEN=454425;CALLER=cnvkit;CORR_CN=8.67;BAF=0.0
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(vcf_content)
            vcf_path = f.name

        try:
            cnvs = get_cnvs(vcf_path)
            all_cnvs = [cnv for callers in cnvs.values() for cnv_list in callers.values() for cnv in cnv_list]
            by_start = {cnv.start: cnv for cnv in all_cnvs}

            # Real, strongly-skewed (negative) measurement -> clamped to a real 0.0.
            assert by_start[10077950].baf == 0.0
            # Raw exact 0.0 -> treated as missing, not a real measurement.
            assert by_start[43475717].baf is None
        finally:
            os.remove(vcf_path)

    def test_evaluate_filter_condition_leaf_operators(self):
        cnv = CNV(
            caller="cnvkit",
            chromosome="chr1",
            genes=["gene1"],
            start=100,
            length=1000,
            type="DUP",
            cn=8.0,
            baf=0.45,
        )

        assert evaluate_filter_condition(cnv, {"field": "adjusted_cn", "operator": ">=", "value": 6.0})
        assert not evaluate_filter_condition(cnv, {"field": "adjusted_cn", "operator": "<", "value": 6.0})
        assert evaluate_filter_condition(cnv, {"field": "baf", "operator": "<", "value": 0.5})
        assert evaluate_filter_condition(cnv, {"field": "caller", "operator": "=", "value": "cnvkit"})
        assert not evaluate_filter_condition(cnv, {"field": "caller", "operator": "=", "value": "Jumble"})
        assert evaluate_filter_condition(cnv, {"field": "caller", "operator": "!=", "value": "Jumble"})
        # A field with no value (e.g. an extra field never extracted for this
        # pipeline, or a genuinely missing length) is always False.
        cnv_missing_length = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=None,
            type="DUP", cn=8.0, baf=0.45,
        )
        assert not evaluate_filter_condition(cnv_missing_length, {"field": "length", "operator": ">", "value": 0})

    def test_evaluate_filter_condition_combinators(self):
        cnv = CNV(
            caller="Jumble", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DEL", cn=1.9, baf=0.0,
        )
        condition = {
            "all_of": [
                {"field": "adjusted_cn", "operator": ">", "value": 1.4},
                {
                    "any_of": [
                        {"all_of": [
                            {"field": "baf", "operator": ">", "value": 0.3},
                            {"field": "baf", "operator": "<", "value": 0.7},
                        ]},
                        {"field": "caller", "operator": "=", "value": "Jumble"},
                    ]
                },
            ]
        }
        # Here baf=0.0 is a real, definite (non-missing) value that fails the
        # range check outright, so this exercises normal two-valued any_of -
        # the caller="Jumble" branch is what makes it true, distinct from the
        # three-valued "unknown propagates" behavior covered below.
        assert evaluate_filter_condition(cnv, condition)
        assert evaluate_filter_condition(cnv, {"not": {"field": "caller", "operator": "=", "value": "cnvkit"}})
        assert not evaluate_filter_condition(cnv, {"not": {"field": "caller", "operator": "=", "value": "Jumble"}})

    def test_evaluate_filter_condition_missing_field_is_unknown_not_false(self):
        # A field with no data (e.g. a caller that never measured BAF for this
        # segment - baf=None, not 0.0) is "unknown", not False. This matters
        # for `not`: not(unknown) must stay unknown rather than flip to True.
        # A naive two-valued evaluator would otherwise say "not BAF-in-neutral-
        # range" is satisfied just because BAF wasn't measured at all - for
        # *any* caller with missing BAF, not only Jumble.
        cnv = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DEL", cn=1.9, baf=None,
        )
        baf_in_neutral_range = {
            "all_of": [
                {"field": "baf", "operator": ">", "value": 0.3},
                {"field": "baf", "operator": "<", "value": 0.7},
            ]
        }
        assert evaluate_filter_condition(cnv, baf_in_neutral_range) is None
        assert evaluate_filter_condition(cnv, {"not": baf_in_neutral_range}) is None

        # "Unknown" mixed into any_of/all_of alongside definite conditions
        # still doesn't produce a definite True classification.
        config = {
            "loh": {
                "all_of": [
                    {"field": "adjusted_cn", "operator": ">=", "value": 1.4},
                    {"field": "adjusted_cn", "operator": "<=", "value": 2.5},
                    {"not": baf_in_neutral_range},
                ]
            }
        }
        assert classify_cnv(cnv, config) is None
        assert not passes_table_filter(cnv, config)

    def test_evaluate_filter_condition_default_substitutes_for_missing_field(self):
        # A leaf may set `default` to opt out of the "missing field is
        # unknown" behavior above, substituting a value instead - only
        # appropriate when a missing value has an obvious neutral stand-in
        # (e.g. no co-located small variant means no elevated-AF signal to
        # worry about, so a missing Twist_AF should read the same as 0).
        cnv = CNV(
            caller="Jumble", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="AMP", cn=8.0, baf=None, extra={},
        )
        twist_af_leaf = {"field": "Twist_AF", "operator": "<=", "value": 0.15, "default": 0}

        # Missing entirely (never extracted for this pipeline, or no
        # co-located variant): default kicks in, leaf resolves definitively.
        assert evaluate_filter_condition(cnv, twist_af_leaf) is True

        # Without `default`, the exact same missing field is unknown, not
        # False - this is the pre-existing behavior `default` opts out of.
        assert evaluate_filter_condition(
            cnv, {"field": "Twist_AF", "operator": "<=", "value": 0.15}
        ) is None

        # A group combining a real value with a missing-but-defaulted field
        # still resolves to a definite match, instead of getting stuck
        # unresolved the way it would without `default` (see the docstring's
        # "amplification" example in docs/cnv_report.md).
        amplification = {
            "all_of": [
                {"field": "adjusted_cn", "operator": ">=", "value": 6.0},
                twist_af_leaf,
            ]
        }
        assert classify_cnv(cnv, {"amplification": amplification}) == "amplification"

        # `default` is only a fallback for missing data - a present value
        # still fails the check normally.
        cnv_with_high_af = CNV(
            caller="Jumble", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="AMP", cn=8.0, baf=None, extra={"Twist_AF": 0.9},
        )
        assert not evaluate_filter_condition(cnv_with_high_af, twist_af_leaf)

    def test_passes_table_filter_group_selection_and_opt_in(self):
        amp_cnv = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DUP", cn=8.0, baf=0.5,
        )
        loss_cnv = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DEL", cn=1.0, baf=0.5,
        )
        config = {
            "amplification": {"all_of": [{"field": "adjusted_cn", "operator": ">=", "value": 6.0}]},
            "loss": {"all_of": [{"field": "adjusted_cn", "operator": "<=", "value": 1.4}]},
        }
        assert passes_table_filter(amp_cnv, config)
        assert passes_table_filter(loss_cnv, config)
        # No config at all -> nothing to filter by, never "passes".
        assert not passes_table_filter(amp_cnv, {})
        assert not passes_table_filter(amp_cnv, None)

    def test_passes_table_filter_generic_extra_field(self):
        # Field names are entirely config-driven — not tied to any specific
        # pipeline's naming (e.g. twist_solid's "Twist_AF"). Here a made-up
        # "cohort_frequency" field demonstrates the mechanism is general.
        cnv = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DUP", cn=8.0, baf=0.5, extra={"cohort_frequency": 0.02},
        )
        config = {
            "amplification": {
                "all_of": [
                    {"field": "adjusted_cn", "operator": ">=", "value": 6.0},
                    {"field": "cohort_frequency", "operator": "<=", "value": 0.15},
                ]
            }
        }
        assert passes_table_filter(cnv, config)

        cnv_common_artifact = CNV(
            caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000,
            type="DUP", cn=8.0, baf=0.5, extra={"cohort_frequency": 0.5},
        )
        assert not passes_table_filter(cnv_common_artifact, config)

    def test_passes_table_filter_real_world_loh_criteria(self):
        # Mirrors twist_solid's real LOH hard-filter criteria (translated to
        # the structured format), and the exact MCPH1 bug-report scenario:
        # all three callers should be excluded. No caller-specific escape
        # clause is needed here - Jumble's BAF is genuinely missing (None,
        # post safe_baf fix) rather than a comparable 0.0, and three-valued
        # logic keeps "unknown BAF" from spuriously satisfying not(...).
        config = {
            "loss": {
                "all_of": [
                    {"not": {"all_of": [
                        {"field": "artifact_af", "operator": ">", "value": 0.15},
                        {"field": "length", "operator": "<", "value": 10_000_000},
                    ]}},
                    {"not": {"any_of": [
                        {"all_of": [
                            {"field": "adjusted_cn", "operator": ">", "value": 1.4},
                            {"all_of": [
                                {"field": "baf", "operator": ">", "value": 0.3},
                                {"field": "baf", "operator": "<", "value": 0.7},
                            ]},
                        ]},
                        {"field": "adjusted_cn", "operator": ">", "value": 2.5},
                    ]}},
                ]
            }
        }

        for caller, cn, baf in [("jumble", 1.97, None), ("cnvkit", 1.85, 0.346), ("gatk", 1.87, 0.448)]:
            cnv = CNV(
                caller=caller, chromosome="chr1", genes=["MCPH1"], start=100, length=1000,
                type="DEL", cn=cn, baf=baf, extra={"artifact_af": 0.01},
            )
            assert not passes_table_filter(cnv, config), f"{caller} should NOT pass the loss filter"

    def test_classify_cnv_arbitrary_named_groups_in_order(self):
        # table_filter_config isn't limited to exactly "amplification"/"loss" —
        # any number of named groups are supported, matched in definition
        # order (first match wins).
        config = {
            "amplification": {"all_of": [{"field": "adjusted_cn", "operator": ">=", "value": 6.0}]},
            "duplication": {"all_of": [
                {"field": "adjusted_cn", "operator": ">", "value": 2.5},
                {"field": "adjusted_cn", "operator": "<", "value": 6.0},
            ]},
            "deletion": {"all_of": [{"field": "adjusted_cn", "operator": "<=", "value": 1.4}]},
        }

        def cnv_with_cn(cn):
            return CNV(caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000, type="", cn=cn, baf=0.5)

        assert classify_cnv(cnv_with_cn(8.0), config) == "amplification"
        assert classify_cnv(cnv_with_cn(5.0), config) == "duplication"
        assert classify_cnv(cnv_with_cn(1.0), config) == "deletion"
        assert classify_cnv(cnv_with_cn(2.0), config) is None  # no group matches -> "copy neutral"

    def test_classify_cnv_filter_false_excludes_group_from_filter_only(self):
        # A group can be marked filter: false to only affect Type
        # classification, without making the row visible under the filter
        # toggle (passes_table_filter uses filter_only=True internally).
        config = {
            "amplification": {"all_of": [{"field": "adjusted_cn", "operator": ">=", "value": 6.0}]},
            "duplication": {
                "filter": False,
                "all_of": [
                    {"field": "adjusted_cn", "operator": ">", "value": 2.5},
                    {"field": "adjusted_cn", "operator": "<", "value": 6.0},
                ],
            },
        }
        cnv = CNV(caller="cnvkit", chromosome="chr1", genes=["gene1"], start=100, length=1000, type="", cn=5.0, baf=0.5)

        # Type classification (filter_only=False) still finds "duplication".
        assert classify_cnv(cnv, config, filter_only=False) == "duplication"
        # But it must not count toward the filter toggle.
        assert classify_cnv(cnv, config, filter_only=True) is None
        assert not passes_table_filter(cnv, config)

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
                                    cn=6,
                                    baf=0.3,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    cn=2,
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
                                    cn=2,
                                    baf=0.5,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene1"],
                                    start=200,
                                    length=100,
                                    type="COPY_NORMAL",
                                    cn=1,
                                    baf=0.5,
                                ),
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene2"],
                                    start=600,
                                    length=100,
                                    type="DEL",
                                    cn=1,
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
                                    cn=6,
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
                                    cn=1,
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
            ),
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
                                    genes=["gene2"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    cn=2,
                                    baf=1.3,
                                ),
                            ],
                        },
                    },
                    # File with deletions
                    {
                        "chr1": {
                            "cnvkit": [
                                # This variant is present in the other file too,
                                # but with annotated with a different gene.
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene3"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    cn=2,
                                    baf=1.3,
                                ),
                            ],
                        },
                    },
                ],
                filtered_cnvs=[
                    # Filtered amplifications
                    defaultdict(lambda: defaultdict(list)),
                    # Filtered deletions
                    {
                        "chr1": {
                            "cnvkit": [
                                CNV(
                                    caller="cnvkit",
                                    chromosome="chr1",
                                    genes=["gene3"],
                                    start=500,
                                    length=100,
                                    type="COPY_NORMAL",
                                    cn=2,
                                    baf=1.3,
                                ),
                            ],
                        },
                    },
                ],
                expected_cnvs={
                    "chr1": {
                        "cnvkit": {
                            "all_count": 1,
                            "pass_filter_count": 1,
                            "genes": [["gene2", "gene3"]],
                        }
                    }
                },
            ),
        ]

        for case in testcases:
            d = merge_cnv_dicts(case.callers, [], [], [], case.chromosomes, case.filtered_cnvs, case.unfiltered_cnvs)

            for chromosome in case.expected_cnvs:
                chromosomedata = [x for x in d if x["chromosome"] == chromosome]
                assert len(chromosomedata) == 1

                for caller in case.expected_cnvs[chromosome]:
                    callerdata = [x for x in chromosomedata[0]["callers"] if x["name"] == caller]
                    assert len(callerdata) == 1

                    exp = case.expected_cnvs[chromosome][caller]

                    assert exp["all_count"] == len(callerdata[0]["cnvs"])
                    for i, c in enumerate(callerdata[0]["cnvs"]):
                        if "genes" not in exp:
                            continue
                        assert len(c.genes) == len(exp["genes"][i])
                        assert all([exp["genes"][i] for gene in c.genes])

                    n_pass_filter = sum(x.passed_filter for x in callerdata[0]["cnvs"])
                    assert n_pass_filter == exp["pass_filter_count"]

    def test_merge_cnv_calls(self):
        unfiltered_files = [
            os.path.join(
                TEST_DIR,
                "../integration/cnv_sv/svdb_query/sample1_T.pathology.svdb_query.annotate_cnv.cnv_amp_genes.vcf",
            ),
            os.path.join(
                TEST_DIR,
                "../integration/cnv_sv/svdb_query/sample1_T.pathology.svdb_query.annotate_cnv.cnv_loh_genes.vcf",
            ),
        ]
        filtered_files = [
            os.path.join(
                TEST_DIR,
                "../integration/cnv_sv/svdb_query/sample1_T.pathology.svdb_query.annotate_cnv.cnv_amp_genes.filter.cnv_hard_filter_amp.vcf",
            ),
            os.path.join(
                TEST_DIR,
                "../integration/cnv_sv/svdb_query/sample1_T.pathology.svdb_query.annotate_cnv.cnv_loh_genes.filter.cnv_hard_filter_loh.vcf",
            ),
        ]

        uf_cnv = []
        f_cnv = []
        for uf_file, f_file in zip(unfiltered_files, filtered_files):
            uf_cnv.append(get_cnvs(uf_file))
            f_cnv.append(get_cnvs(f_file))

        merged = merge_cnv_calls(uf_cnv, f_cnv)

        cnvkit_cnvs = [x for x in merged if x.caller == "cnvkit"]
        gatk_cnvs = [x for x in merged if x.caller == "gatk"]
        jumble_cnvs = [x for x in merged if x.caller == "jumble"]

        assert len(cnvkit_cnvs) == 7 + 1  # amp + loh
        assert len(gatk_cnvs) == 2
        assert len(jumble_cnvs) == 1

        assert sum(x.passed_filter for x in cnvkit_cnvs) == 4 + 1  # 3 amp + 1 del + 1 overlap with gatk
        assert sum(x.passed_filter for x in gatk_cnvs) == 1 + 1  # 1 amp + 1 for overlap with cnvkit
        assert sum(x.passed_filter for x in jumble_cnvs) == 1  # 1 amp

    def test_cnv_sorting(self):
        """
        Make sure that numerical chromosome names are sorted properly.
        """
        cnvs = [
            CNV(caller="cnvkit", chromosome="chr10", genes=[], start=100, length=200, cn=2, baf=0.5, type="COPY_NORMAL"),
            CNV(caller="cnvkit", chromosome="chr2", genes=[], start=500, length=200, cn=2, baf=0.5, type="COPY_NORMAL"),
            CNV(caller="cnvkit", chromosome="chr1", genes=[], start=500, length=200, cn=2, baf=0.5, type="COPY_NORMAL"),
            CNV(caller="cnvkit", chromosome="chr1", genes=[], start=100, length=200, cn=2, baf=0.5, type="COPY_NORMAL"),
            CNV(caller="cnvkit", chromosome="chrM", genes=[], start=100, length=200, cn=2, baf=0.5, type="COPY_NORMAL"),
        ]

        sorted_cnvs = sort_cnvs(cnvs)

        assert len(sorted_cnvs) == len(cnvs)
        assert sorted_cnvs[0] == cnvs[3]
        assert sorted_cnvs[1] == cnvs[2]
        assert sorted_cnvs[2] == cnvs[1]
        assert sorted_cnvs[3] == cnvs[0]
        assert sorted_cnvs[4] == cnvs[4]
