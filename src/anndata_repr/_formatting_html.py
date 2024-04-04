from __future__ import annotations

import typing
import uuid
from html import escape
from anndata_repr._create_icons import get_display

from anndata_repr._formatting_table import dataframe_to_table

from ._formatting_dask_svg import svg_2d
from ._formatting_html_xarray import (
    _icon,
    _obj_repr,
    collapsible_section,
    short_data_repr_html,
    inline_variable_array_repr,
)
import pandas as pd

if typing.TYPE_CHECKING:
    import anndata
    from anndata._core.aligned_mapping import PairwiseArraysBase, AxisArraysBase

__all__ = ["format_anndata_html"]


def summarize_series(obj: pd.Series) -> str:
    """Summarize attributes of Pandas Series.
    List of dtype, number of elements, number of unique elements, number of not None elements.
    """
    len_values = len(obj)
    unique_values = len(obj.unique())
    nonnull_values = obj.count()
    dtype_values = obj.dtype
    enum = {
        "dtype": dtype_values,
        "Items": len_values,
        "Unique items": unique_values,
        "Non-null items": nonnull_values,
    }
    attrs_dl = "".join(
        f"<dt><span>{escape(str(k))} :</span></dt><dd>{escape(str(v))}</dd>"
        for k, v in enum.items()
    )
    return f"<dl class='ad-attrs-data'>{attrs_dl}</dl>"


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
    cssclass_idx = " class='ad-has-index'" if is_index else ""

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())

    preview = escape(inline_variable_array_repr(col, 35))
    attrs_ul = summarize_series(col) if isinstance(col, pd.Series) else ""
    data_repr = short_data_repr_html(col)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    shape = "" if len(col.shape) == 1 else f"({', '.join(map(str, col.shape))})"

    return (
        f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='ad-var-dims'>{shape}</div>"
        f"<div class='ad-var-dtype'>{dtype}</div>"
        f"<div class='ad-var-preview ad-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='ad-var-attrs-in' type='checkbox'>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>{attrs_icon}</label>"
        f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>{data_icon}</label>"
        f"<div class='ad-var-attrs'>{attrs_ul}</div>"
        f"<div class='ad-var-data'>{data_repr}</div>"
    )



def summarize_table(df: pd.DataFrame, is_index: bool = True) -> str:
    name = "Table"
    cssclass_idx = " class='ad-has-index'" if is_index else ""

    # "unique" ids required to expand/collapse subsections
    data_id = "data-" + str(uuid.uuid4())

    data_repr = dataframe_to_table(df[0:10])

    data_icon = _icon("icon-database")

    return (
        f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>{data_icon}</label>"
        f"<div class='ad-var-data'>{data_repr}</div>"
    )


def summarize_obsvar(obsvar: pd.DataFrame) -> str:
    """Summarize the obs or var DataFrame.

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
    for k in obsvar:
        assert isinstance(k, str), "Column of dataframe is not a string"
        li_content = summarize_columns(k, obsvar[k])
        li_items.append(f"<li class='ad-var-item'>{li_content}</li>")
    vars_li = "".join(li_items)
    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def summarize_arrays(arr: PairwiseArraysBase | AxisArraysBase) -> str:
    li_items = []
    for k, v in arr.items():
        li_content = summarize_columns(k, v)
        li_items.append(f"<li class='ad-var-item'>{li_content}</li>")
    vars_li = "".join(li_items)
    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def summarize_layer(name: str, layer, is_index: bool = True) -> str:
    name = escape(str(name))
    dtype = escape(str(layer.dtype))
    cssclass_idx = " class='ad-has-index'" if is_index else ""

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())

    preview = escape(inline_variable_array_repr(layer, 35))
    attrs_ul = ""  # summarize_attrs(col)
    data_repr = short_data_repr_html(layer)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='ad-var-dims'>(X)</div>"
        f"<div class='ad-var-dtype'>{dtype}</div>"
        f"<div class='ad-var-preview ad-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='ad-var-attrs-in' type='checkbox'>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>{attrs_icon}</label>"
        f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>{data_icon}</label>"
        f"<div class='ad-var-attrs'>{attrs_ul}</div>"
        f"<div class='ad-var-data'>{data_repr}</div>"
    )


def summarize_X(adata: anndata.AnnData) -> str:
    li_items = []
    for layer_name in adata.layers.keys():
        li_content = summarize_layer(layer_name, adata.layers[layer_name])
        li_items.append(f"<li class='ad-var-item'>{li_content}</li>")
    vars_li = "".join(li_items)
    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def array_section(adata) -> str:
    display = get_display(adata)
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())
    collapsed = True
    # TODO: Always use the svg_2d for the preview?
    def convert_newlines_to_br(html_string):
        return html_string.replace('\n', '<br>')
    data_repr = f"<p>{convert_newlines_to_br(adata.__repr__())}</p>"#short_data_repr_html(X)
    data_icon = _icon("icon-database")


    return (
        "<div class='ad-array-wrap version3'>"
        f"<input id='{data_id}' class='ad-array-in' type='checkbox' {collapsed}>"
        f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
        f"<div class='ad-array-preview ad-preview'><span>{display}</span></div>"
        f"<div class='ad-array-data'>{data_repr}</div>"
        "</div>"
    )


def summaize_uns(uns: dict) -> str:
    li_items = []

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    def summarize_un(name, un):
        is_index = True
        cssclass_idx = " class='ad-has-index'" if is_index else ""
        # "unique" ids required to expand/collapse subsections
        attrs_id = "attrs-" + str(uuid.uuid4())
        data_id = "data-" + str(uuid.uuid4())
        attrs_ul = ""
        data_repr = f"<pre>{escape(repr(un))}</pre>"
        type_name = type(un).__name__
        return (
            f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
            f"<div class='ad-var-dims'></div>"
            f"<div class='ad-var-dtype'></div>"
            f"<div class='ad-var-preview ad-preview'>{type_name}</div>"
            f"<input id='{attrs_id}' class='ad-var-attrs-in' type='checkbox'>"
            f"<label for='{attrs_id}' title='Show/Hide attributes'>{attrs_icon}</label>"
            f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
            f"<label for='{data_id}' title='Show/Hide data repr'>{data_icon}</label>"
            f"<div class='ad-var-attrs'>{attrs_ul}</div>"
            f"<div class='ad-var-data'>{data_repr}</div>"
        )

    for k, v in uns.items():
        li_content = summarize_un(k, v)
        li_items.append(f"<li class='ad-var-item'>{li_content}</li>")
    vars_li = "".join(li_items)
    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def format_anndata_html(adata: anndata.AnnData) -> str:
    """Format an AnnData object as an HTML string.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to format.

    """
    obj_type = f"anndata.{type(adata).__name__}"
    arr_name = ""  # TODO: add somethign here?

    dims_li = "".join(
        f"<li><span class='ad-has-index'>obs</span>: {adata.n_obs}</li>"
        f"<li><span class='ad-has-index'>var</span>: {adata.n_vars}</li>"
    )

    header_components = [
        f"<div class='ad-obj-type'>{obj_type}</div>",
        f"<div class='ad-array-name'>{arr_name}</div>",
        f"<ul class='ad-dim-list'>{dims_li}</ul>",
    ]


    sections = [
        array_section(adata),
        collapsible_section(
            "layers",
            details=summarize_X(adata),
            n_items=len(adata.layers),
            collapsed=True,
        ),
        collapsible_section(
            "obs",
            details=summarize_obsvar(adata.obs),
            n_items=len(adata.obs.columns),
            collapsed=True,
        ),
        collapsible_section(
            "var",
            details=summarize_obsvar(adata.var),
            n_items=len(adata.var.columns),
            collapsed=True,
        ),
        collapsible_section(
            "obsm",
            details=summarize_arrays(adata.obsm),
            n_items=len(adata.obsm),
            collapsed=True,
        ),
        collapsible_section(
            "obsp",
            details=summarize_arrays(adata.obsp),
            n_items=len(adata.obsp),
            collapsed=True,
        ),
        collapsible_section(
            "varm",
            details=summarize_arrays(adata.varm),
            n_items=len(adata.varm),
            collapsed=True,
        ),
        collapsible_section(
            "varp",
            details=summarize_arrays(adata.varp),
            n_items=len(adata.varp),
            collapsed=True,
        ),
        collapsible_section(
            "uns",
            details=summaize_uns(adata.uns),
            n_items=len(adata.uns),
            collapsed=True,
        ),
    ]

    return _obj_repr(adata, header_components, sections)
