import pathlib
import pytest
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from workflow.scripts.general_html_report import *


def test_valid_table():
    table = [
        {
            "col1": "val1",
            "col2": 1,
        },
        {
            "col1": "val2",
            "col2": 2,
        },
        {
            "col1": "val3",
            "col2": 3,
        },
    ]

    assert validate_table_data(table)


def test_empty_table():
    with pytest.raises(ValueError, match="empty table"):
        validate_table_data([])


def test_unequal_columns():
    table = [
        {
            "col1": "val1",
            "col2": 1,
        },
        {
            "col1": "val2",
            "col3": 2,
        },
        {
            "col1": "val3",
            "col2": 3,
        },
    ]

    with pytest.raises(
        ValueError,
        match="expected columns 'col1', 'col2' in row 2, found 'col1', 'col3'",
    ):
        validate_table_data(table)


def test_missing_columns():
    table = [
        {
            "col1": "val1",
            "col2": 1,
        },
        {
            "col1": "val2",
        },
        {
            "col1": "val3",
            "col2": 3,
        },
    ]

    with pytest.raises(ValueError, match="expected 2 columns in row 2, found 1"):
        validate_table_data(table)
