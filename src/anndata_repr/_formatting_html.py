from __future__ import annotations

import math
import typing
import uuid
from html import escape
from itertools import zip_longest

import numpy as np
from xarray.core.formatting import first_n_items, format_items, last_n_items

from ._formatting_html_xarray import (
    _icon,
    _obj_repr,
    collapsible_section,
    short_data_repr_html,
)

if typing.TYPE_CHECKING:
    import anndata
    import pandas as pd

__all__ = ["format_anndata_html"]


# from xarray.core.formatting import format_array_flat
def format_array_flat(array: pd.Series, max_width: int):
    """Return a formatted string for as many items in the flattened version of
    array that will fit within max_width characters.
    """
    # every item will take up at least two characters, but we always want to
    # print at least first and last items
    max_possibly_relevant = min(max(array.size, 1), max(math.ceil(max_width / 2.0), 2))
    relevant_front_items = format_items(
        first_n_items(array, (max_possibly_relevant + 1) // 2)
    )
    relevant_back_items = format_items(last_n_items(array, max_possibly_relevant // 2))
    # interleave relevant front and back items:
    #     [a, b, c] and [y, z] -> [a, z, b, y, c]
    relevant_items = sum(
        zip_longest(relevant_front_items, reversed(relevant_back_items)), ()
    )[:max_possibly_relevant]

    cum_len = np.cumsum([len(s) + 1 for s in relevant_items]) - 1
    if (array.size > 2) and (
        (max_possibly_relevant < array.size) or (cum_len > max_width).any()
    ):
        padding = " ... "
        max_len = max(int(np.argmax(cum_len + len(padding) - 1 > max_width)), 2)
        count = min(array.size, max_len)
    else:
        count = array.size
        padding = "" if (count <= 1) else " "

    num_front = (count + 1) // 2
    num_back = count - num_front
    # note that num_back is 0 <--> array.size is 0 or 1
    #                         <--> relevant_back_items is []
    pprint_str = "".join(
        [
            " ".join(relevant_front_items[:num_front]),
            padding,
            " ".join(relevant_back_items[-num_back:]),
        ]
    )

    # As a final check, if it's still too long even with the limit in values,
    # replace the end with an ellipsis
    # NB: this will still returns a full 3-character ellipsis when max_width < 3
    if len(pprint_str) > max_width:
        pprint_str = pprint_str[: max(max_width - 3, 0)] + "..."

    return pprint_str


def maybe_truncate(obj, maxlen=500):
    s = str(obj)
    if len(s) > maxlen:
        s = s[: (maxlen - 3)] + "..."
    return s


def inline_variable_array_repr(col: pd.Series, max_width: int):
    """Build a one-line summary of a variable's data."""
    return format_array_flat(col, max_width)


def summarize_columns(name: str, col: pd.Series) -> str:
    """Summarize a single column of a DataFrame.

    Parameters
    ----------
    name : str
        The name of the column.

    col : pd.Series
        The column to summarize.

    Returns
    -------
    str
        The HTML representation of the column.
    """
    name = escape(str(name))
    dtype = escape(str(col.dtype))

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())
    # disabled = "" if len(var.attrs) else "disabled"

    preview = escape(inline_variable_array_repr(col, 35))
    # attrs_ul = summarize_attrs(var.attrs)
    data_repr = short_data_repr_html(col)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='xr-var-dtype'>{dtype}</div>"
        f"<div class='xr-var-preview xr-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='xr-var-attrs-in' "
        f"type='checkbox'>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>"
        f"{attrs_icon}</label>"
        f"<input id='{data_id}' class='xr-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>"
        f"{data_icon}</label>"
        # f"<div class='xr-var-attrs'>{attrs_ul}</div>"
        f"<div class='xr-var-data'>{data_repr}</div>"
    )


def summarize_obs(variables: pd.DataFrame) -> str:
    """Summarize the variables of a DataFrame.

    Parameters
    ----------
    variables : pd.DataFrame
        The obs or var DataFrame to summarize.

    Returns
    -------
    str
        The HTML representation of the variables.
    """
    li_items = []
    for k in variables:
        li_content = summarize_columns(k, typing.cast("pd.Series", variables[k]))
        li_items.append(f"<li class='xr-var-item'>{li_content}</li>")

    vars_li = "".join(li_items)

    return f"<ul class='xr-var-list'>{vars_li}</ul>"


def format_anndata_html(adata: anndata.AnnData) -> str:
    """Format an AnnData object as an HTML string.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to format.

    """
    obj_type = f"{type(adata).__name__}"

    header_components = [f"<div class='xr-obj-type'>{escape(obj_type)}</div>"]

    sections = [
        collapsible_section(
            "obs",
            details=summarize_obs(adata.obs),
            n_items=3,
            enabled=True,
            collapsed=True,
        ),
        collapsible_section(
            "var",
            details=summarize_obs(adata.var),
            n_items=3,
            enabled=True,
            collapsed=True,
        ),
    ]

    return _obj_repr(adata, header_components, sections)
