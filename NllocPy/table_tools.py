# =============================================================================
# Author: Luigi Sante Zampa
# Affiliation: NllocPy project (institutional affiliation not specified)
# Created: 2026-11-04
# Last modified: 2026-11-04
#
# Purpose:
#   Convert tables into dictionaries of aligned, one-dimensional NumPy arrays.
#   All computational code can therefore use dictionary-based data structures,
#   while notebooks may supply/receive pandas DataFrames when requested.
#   Pandas remains optional; it is not an installation requirement of NllocPy.
# =============================================================================

# -----------------------------------------------------------------------------
# Import external libraries
# -----------------------------------------------------------------------------
from collections.abc import Mapping

import numpy as np

try:
    import pandas as pd
except ImportError:  # Keep NllocPy usable in environments without pandas.
    pd = None


# -----------------------------------------------------------------------------
# Shared module constants
# -----------------------------------------------------------------------------
_DATAFRAME_UNAVAILABLE = "pandas is required to return a DataFrame"


# =============================================================================
# Input format detection and conversion
# =============================================================================

# -----------------------------------------------------------------------------
def is_dataframe(data):
    """Identify a pandas DataFrame without requiring pandas to be installed.

    Parameters
    ----------
    data : object
        Any Python object.

    Returns
    -------
    bool
        True only when pandas is available and the object is a DataFrame.

    Examples
    --------
    >>> is_dataframe({"x": [1, 2]})
    False
    """
    output = pd is not None and isinstance(data, pd.DataFrame)
    return output


# -----------------------------------------------------------------------------
def to_field_dict(data, *, copy=True):
    """Convert a mapping or DataFrame to aligned one-dimensional NumPy fields.

    DataFrames are converted immediately. Every dictionary field is retained;
    the function does not silently reorder samples or change their values.

    Parameters
    ----------
    data : dict or pandas.DataFrame
        Table with a field/column for each one-dimensional array.
    copy : bool, default True
        If True, copy the arrays so later modifications cannot change input.

    Returns
    -------
    dict
        Dictionary mapping field names to equally sized NumPy arrays.

    Raises
    ------
    TypeError
        Input is not a mapping or DataFrame.
    ValueError
        A column is not one-dimensional or columns have different lengths.

    Examples
    --------
    >>> fields = to_field_dict({"x": [2, 1], "vp": [5.0, 4.0]})
    >>> fields["vp"].tolist()
    [5.0, 4.0]
    """
    if is_dataframe(data):
        original = data.to_dict(orient="list")
    elif isinstance(data, Mapping):
        original = dict(data)
    else:
        raise TypeError("Expected a dictionary-like mapping or pandas DataFrame")

    fields = {}
    expected_length = None

    # Explicitly validate columns so they remain correctly aligned.
    for key, values in original.items():
        array = np.asarray(values)

        if array.ndim != 1:
            raise ValueError(f"Field {key!r} must be one-dimensional")

        if expected_length is None:
            expected_length = array.size

        if array.size != expected_length:
            raise ValueError("All fields must contain the same number of rows")

        fields[key] = array.copy() if copy else array

    output = fields
    return output


# =============================================================================
# Optional reconstruction of the caller's table format
# =============================================================================

# -----------------------------------------------------------------------------
def restore_field_table(fields, *, template=None, as_dataframe=False):
    """Return a dictionary, or reconstruct a pandas DataFrame when requested.

    Use this helper only at the interface boundary of a function. The
    numerical processing between input and output should use dictionaries.

    Parameters
    ----------
    fields : dict
        One-dimensional, equally sized output fields.
    template : pandas.DataFrame or None, optional
        When a DataFrame is passed, return a DataFrame. Its index is reused
        only when its number of rows equals that of the output fields.
    as_dataframe : bool, default False
        Force a pandas DataFrame even without a DataFrame template.

    Returns
    -------
    dict or pandas.DataFrame
        Dictionary by default; DataFrame if requested or originally supplied.

    Raises
    ------
    ImportError
        A DataFrame is requested but pandas is not installed.

    Examples
    --------
    >>> restore_field_table({"x": np.array([1.0])})["x"].tolist()
    [1.0]
    """
    converted = to_field_dict(fields)
    wants_frame = as_dataframe or is_dataframe(template)

    if not wants_frame:
        output = converted
        return output

    if pd is None:
        raise ImportError(_DATAFRAME_UNAVAILABLE)

    length = len(next(iter(converted.values()))) if converted else 0
    index = None

    # Preserve the caller's index only if its length remains compatible.
    if is_dataframe(template) and len(template.index) == length:
        index = template.index.copy()

    output = pd.DataFrame(converted, index=index)
    return output
