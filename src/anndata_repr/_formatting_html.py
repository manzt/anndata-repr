from __future__ import annotations

import typing
import uuid
from html import escape


from ._svg import anndata_svg, _COLORS
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


def dataframe_to_table(dataframe, max_len=10):
    # Initialize the table with the correct class names for styling
    table = """<div class='relative overflow-x-auto my-div'>
    <table class="">
    """

    def make_header(columns):
        # Style the header row according to the provided CSS class names
        header = '<thead class="table-header"><tr>'
        for column in columns:
            header += '<th scope="col" class="column-header">' + column + "</th>"
        header += "</tr></thead>"
        return header

    def make_truncated_data(data, max_len):
        # Function to truncate data if it's longer than 10 characters
        if len(str(data)) > max_len:
            return str(data)[:max_len] + "..."
        return str(data)

    def make_row(row, index, max_len):
        # Style each row according to the provided CSS class names, alternating row color not implemented in CSS
        row_html = f'<tr class="table-row">'
        for value in row:
            row_html += (
                f'<td class="table-cell">{make_truncated_data(value, max_len)}</td>'
            )
        row_html += "</tr>"
        return row_html

    # Construct the table header
    table += make_header(dataframe.columns)

    # Construct each row of the table
    for index, row in enumerate(dataframe.itertuples(index=False), start=1):
        table += make_row(row, index, max_len)

    # Close the table and div tags
    table += "</table></div>"
    return table


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


def summarize_obsvar(obsvar: pd.DataFrame, as_df: bool = True) -> str:
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
    if as_df:
        return obsvar._repr_html_()

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


def array_section(adata: anndata.AnnData) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())
    collapsed = True
    preview = anndata_svg(adata)
    data_repr = f'<pre>{escape(repr(adata))}</pre>'
    data_icon = _icon("icon-database")

    return (
        "<div class='ad-array-wrap'>"
        f"<input id='{data_id}' class='ad-array-in' type='checkbox' {collapsed}>"
        f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
        f"<div class='ad-array-preview ad-preview'><span>{preview}</span></div>"
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
            color=_COLORS["X"].hover,
            collapsed=True,
        ),
        collapsible_section(
            "obs",
            details=summarize_obsvar(adata.obs, as_df=True),
            n_items=len(adata.obs.columns),
            color=_COLORS["obs"].hover,
            collapsed=True,
        ),
        collapsible_section(
            "var",
            details=summarize_obsvar(adata.var, as_df=True),
            n_items=len(adata.var.columns),
            color=_COLORS["var"].hover,
            collapsed=True,
        ),
        collapsible_section(
            "obsm",
            details=summarize_arrays(adata.obsm),
            n_items=len(adata.obsm),
            color=_COLORS["obsm"].hover,
            collapsed=True,
        )
        if len(adata.obsm)
        else "",
        collapsible_section(
            "obsp",
            details=summarize_arrays(adata.obsp),
            n_items=len(adata.obsp),
            color=_COLORS["obsp"].hover,
            collapsed=True,
        )
        if len(adata.obsp)
        else "",
        collapsible_section(
            "varm",
            details=summarize_arrays(adata.varm),
            n_items=len(adata.varm),
            color=_COLORS["varm"].hover,
            collapsed=True,
        )
        if len(adata.varm)
        else "",
        collapsible_section(
            "varp",
            details=summarize_arrays(adata.varp),
            n_items=len(adata.varp),
            color=_COLORS["varp"].hover,
            collapsed=True,
        )
        if len(adata.varp)
        else "",
        collapsible_section(
            "uns",
            details=summaize_uns(adata.uns),
            n_items=len(adata.uns),
            collapsed=True,
        )
        if len(adata.uns)
        else "",
    ]

    return _obj_repr(adata, header_components, sections)
