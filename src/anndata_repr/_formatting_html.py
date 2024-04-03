from __future__ import annotations

import math
import typing
import uuid
from html import escape
from itertools import zip_longest

import numpy as np
import pandas as pd
from xarray.core.formatting import (
    first_n_items,
    format_items,
    last_n_items,
    short_array_repr,
)

from ._formatting_html_xarray import (
    _icon,
    _obj_repr,
    collapsible_section,
)

from ._formatting_dask_svg import svg_2d

if typing.TYPE_CHECKING:
    import anndata
    import pandas as pd

__all__ = ["format_anndata_html"]

def summarize_attrs(obj: pd.Series) -> str:
    """Summarize attributes of Pandas Series.
    List of dtype, number of elements, number of unique elements, number of not None elements.
    """
    len_values = len(obj)
    unique_values = len(obj.unique())
    nonnull_values = obj.count()
    dtype_values = obj.dtype
    enum = {'dtype': dtype_values,
             'Items': len_values, 
             'Unique items': unique_values,
             'Non-null items': nonnull_values}
             
    attrs_dl = "".join(
        f"<dt><span>{escape(str(k))} :</span></dt><dd>{escape(str(v))}</dd>"
        for k, v in enum.items()
    )
    return f"<dl class='xr-attrs-data'>{attrs_dl}</dl>"


def short_data_repr_html(obj: pd.Series) -> str:
    if hasattr(obj, "_repr_html_"):
        return obj._repr_html_()
    text = short_array_repr(obj)
    return f"<pre>{escape(text)}</pre>"


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
            " ".join(relevant_front_items[:num_front]),  # type: ignore
            padding,
            " ".join(relevant_back_items[-num_back:]),  # type: ignore
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


def summarize_category(col):
    counts = col.value_counts()
    counts = counts.reset_index()
    counts.columns = ["value", "count"]
    counts = counts.sort_values(by="count", ascending=False)

    if len(counts) > 10:
        counts = counts.head(10) # take top 10
        # counts = pd.concat([counts,(
        #     pd.DataFrame({"value": "...", "count": len(col) - counts["count"].sum()})
        # )])

    counts_html = "".join(
        f"<dt>{escape(str(row['value']))}</dt><dd>{row['count']}</dd>"
        for _, row in counts.iterrows()
    )
    return counts_html

# def histogram_to_svg(hist):



def summarize_numeric(col):
    # get min and max of pandas series, compute bins and draw and svg histogram
    min_val = col.min()
    max_val = col.max()
    bins = np.linspace(min_val, max_val, 10)
    hist, _ = np.histogram(col, bins=bins)
    hist = hist / hist.sum()
    print(hist)
    # svg_content = histogram_to_svg(hist)

    return hist


def summarize_value_counts(name: str, col: pd.Series, is_index: bool = True) -> str:
    """Summarize the value counts of a Series."""

    # if column is category, show value counts 
    if col.dtype.name == "category":
        counts_html = summarize_category(col)
    elif col.dtype.name == "object":
        counts_html = summarize_category(col)
    elif col.dtype.name == "bool":
        counts_html = summarize_category(col)
    elif col.dtype.name == "float32":
        counts_html = summarize_numeric(col)
    elif col.dtype.name == "int64":
        counts_html = summarize_numeric(col)


        
    

    return f"<dl class='xr-value-counts'>{counts_html}</dl>"




def summarize_columns(name: str, col: pd.Series, is_index: bool = True) -> str:
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
    cssclass_idx = " class='xr-has-index'" if is_index else ""

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())

    preview = escape(inline_variable_array_repr(col, 35))
    attrs_ul = summarize_attrs(col)
    data_repr = short_data_repr_html(col)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    data_summary = summarize_value_counts(name, col, is_index)

    return (
        f"<div class='xr-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='xr-var-dims'></div>"
        f"<div class='xr-var-dtype'>{dtype}</div>"
        f"<div class='xr-var-preview xr-preview'>{preview}</div>"
        f"<div class='xr-var-summary' styles='display:none'>{data_summary}</div>"
        f"<input id='{attrs_id}' class='xr-var-attrs-in' "
        f"type='checkbox'>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>"
        f"{attrs_icon}</label>"
        f"<input id='{data_id}' class='xr-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>"
        f"{data_icon}</label>"
        f"<div class='xr-var-attrs'>{attrs_ul}</div>"
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
        assert isinstance(k, str), "Column of dataframe is not a string"
        li_content = summarize_columns(k, variables[k])
        li_items.append(f"<li class='xr-var-item'>{li_content}</li>")

    vars_li = "".join(li_items)

    return f"<ul class='xr-var-list'>{vars_li}</ul>"


def format_var_obs(adata: anndata.AnnData) -> str:
    dims_li = "".join(
        f"<li><span class='xr-has-index'>obs</span>: {adata.n_obs}</li>"
        f"<li><span class='xr-has-index'>var</span>: {adata.n_vars}</li>"
    )
    return f"<ul class='xr-dim-list'>{dims_li}</ul>"


def array_section(X) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())
    collapsed = True
    # TODO: Always use the svg_2d for the preview?
    preview = svg_2d((tuple((dim,) for dim in X.shape)))
    data_repr = short_data_repr_html(X)
    data_icon = _icon("icon-database")

    return (
        "<div class='xr-array-wrap'>"
        f"<input id='{data_id}' class='xr-array-in' type='checkbox' {collapsed}>"
        f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
        f"<div class='xr-array-preview xr-preview'><span>{preview}</span></div>"
        f"<div class='xr-array-data'>{data_repr}</div>"
        "</div>"
    )


def format_anndata_html(adata: anndata.AnnData) -> str:
    """Format an AnnData object as an HTML string.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to format.

    """
    obj_type = f"anndata.{type(adata).__name__}"
    arr_name = ""  # TODO: add somethign here?

    header_components = [
        f"<div class='xr-obj-type'>{obj_type}</div>",
        f"<div class='xr-array-name'>{arr_name}</div>",
        format_var_obs(adata),
    ]

    sections = [
        array_section(adata.X),
        collapsible_section(
            "obs",
            details=summarize_obs(adata.obs),
            n_items=len(adata.obs.columns),
            enabled=True,
            collapsed=False,
        ),
        collapsible_section(
            "var",
            details=summarize_obs(adata.var),
            n_items=len(adata.var.columns),
            enabled=True,
            collapsed=True,
        ),
    ]

    return _obj_repr(adata, header_components, sections)
