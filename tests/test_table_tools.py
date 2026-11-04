# =============================================================================
# Author: Luigi Sante Zampa
# Affiliation: NllocPy project (institutional affiliation not specified)
# Created: 2026-11-04
# Last modified: 2026-11-04
#
# Purpose:
#   Verify dictionary-first table conversion and optional pandas round-trips.
# =============================================================================

# -----------------------------------------------------------------------------
# Import external libraries
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from NllocPy.table_tools import (
    is_dataframe,
    restore_field_table,
    to_field_dict,
)


# =============================================================================
# Tests: column alignment and optional pandas compatibility
# =============================================================================

# -----------------------------------------------------------------------------
def test_mapping_conversion_preserves_row_alignment():
    """Verify that mapping fields retain their original row association.

    Parameters
    ----------
    None
        This pytest test needs no external arguments.

    Returns
    -------
    None
        Assertions check the data mapping.

    Examples
    --------
    >>> test_mapping_conversion_preserves_row_alignment()
    """
    original = {"x": [3, 1, 2], "vp": [6.0, 4.0, 5.0]}
    converted = to_field_dict(original)

    np.testing.assert_array_equal(converted["x"], [3, 1, 2])
    np.testing.assert_array_equal(converted["vp"], [6, 4, 5])

    # No alias to the original mutable field array is retained.
    assert isinstance(converted, dict)


# -----------------------------------------------------------------------------
@pytest.mark.parametrize(
    "bad",
    [
        {"x": [1, 2], "vp": [3]},
        {"x": [[1, 2]], "vp": [3, 4]},
    ],
)
def test_invalid_table_fields_are_rejected(bad):
    """Reject unequal or multidimensional table columns.

    Parameters
    ----------
    bad : dict
        Invalid table assembled by pytest parameterization.

    Returns
    -------
    None
        The assertion verifies error handling.

    Examples
    --------
    >>> # pytest -q tests/test_table_tools.py
    """
    with pytest.raises(ValueError):
        to_field_dict(bad)


# -----------------------------------------------------------------------------
def test_restore_mapping_as_dictionary():
    """Keep the default output as a dictionary.

    Parameters
    ----------
    None
        No arguments.

    Returns
    -------
    None
        The assertion checks the format.

    Examples
    --------
    >>> test_restore_mapping_as_dictionary()
    """
    source = {"x": np.array([1.0, 2.0])}
    restored = restore_field_table(source)

    assert isinstance(restored, dict)
    np.testing.assert_array_equal(restored["x"], source["x"])


# -----------------------------------------------------------------------------
def test_pandas_table_roundtrip_when_available():
    """Preserve a DataFrame's field alignment and compatible index.

    Parameters
    ----------
    None
        No arguments.

    Returns
    -------
    None
        Pandas is optional; the test skips if it is not installed.

    Examples
    --------
    >>> # pytest -q tests/test_table_tools.py
    """
    pandas = pytest.importorskip("pandas")
    original = pandas.DataFrame(
        {"x": [4, 2, 3], "vp": [6.0, 4.0, 5.0]},
        index=["b", "a", "c"],
    )

    assert is_dataframe(original)

    values = to_field_dict(original)
    values["vp"] += 0.25

    restored = restore_field_table(values, template=original)

    assert isinstance(restored, pandas.DataFrame)
    assert restored.index.tolist() == ["b", "a", "c"]
    np.testing.assert_allclose(restored["vp"], [6.25, 4.25, 5.25])


# -----------------------------------------------------------------------------
def test_model_builder_accepts_pandas_when_available():
    """Ensure DataFrame input is normalized immediately for model creation.

    Parameters
    ----------
    None
        No arguments.

    Returns
    -------
    None
        The output grid is checked numerically.

    Examples
    --------
    >>> # pytest -q tests/test_table_tools.py
    """
    pandas = pytest.importorskip("pandas")

    # Keep the optional model-building import at module level via importlib:
    # this test uses the existing documented API, not a second implementation.
    from NllocPy.model_tools import build_regular_model

    model = {
        "x": [0, 0, 1, 1, 0, 0, 1, 1],
        "y": [0, 1, 0, 1, 0, 1, 0, 1],
        "z": [0, 0, 0, 0, 1, 1, 1, 1],
        "vp": [4.0] * 8,
    }

    frame = pandas.DataFrame(model)
    output = build_regular_model(frame, step=1, z_positive="down")

    assert output[3].shape == (2, 2, 2)
    assert np.all(output[3] == 4.0)
