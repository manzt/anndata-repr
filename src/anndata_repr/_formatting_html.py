from __future__ import annotations

import typing
import uuid
from html import escape

from ._formatting_dask_svg import svg_2d
from ._formatting_html_xarray import (
    _icon,
    _obj_repr,
    collapsible_section,
    short_data_repr_html,
    inline_variable_array_repr,
)

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

    return (
        f"<div class='xr-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='xr-var-dims'></div>"
        f"<div class='xr-var-dtype'>{dtype}</div>"
        f"<div class='xr-var-preview xr-preview'>{preview}</div>"
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
