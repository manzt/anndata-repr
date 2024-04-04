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
    from anndata._core.aligned_mapping import LayersBase
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
    attrs_ul = summarize_attrs(col)
    data_repr = short_data_repr_html(col)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='ad-var-dims'></div>"
        f"<div class='ad-var-dtype'>{dtype}</div>"
        f"<div class='ad-var-preview ad-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='ad-var-attrs-in' type='checkbox'>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>{attrs_icon}</label>"
        f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>{data_icon}</label>"
        f"<div class='ad-var-attrs'>{attrs_ul}</div>"
        f"<div class='ad-var-data'>{data_repr}</div>"
    )


def dataframe_to_table(dataframe):
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

    def make_truncated_data(data):
        # Function to truncate data if it's longer than 10 characters
        if len(str(data)) > 10:
            return str(data)[:10] + "..."
        return str(data)

    def make_row(row, index):
        # Style each row according to the provided CSS class names, alternating row color not implemented in CSS
        row_html = f'<tr class="table-row">'
        for value in row:
            row_html += f'<td class="table-cell">{make_truncated_data(value)}</td>'
        row_html += "</tr>"
        return row_html

    # Construct the table header
    table += make_header(dataframe.columns)

    # Construct each row of the table
    for index, row in enumerate(dataframe.itertuples(index=False), start=1):
        table += make_row(row, index)

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
        li_items.append(f"<li class='ad-var-item'>{li_content}</li>")

    vars_li = "".join(li_items)

    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def summarize_layer(name: str, layer, is_index: bool = True) -> str:
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


def array_section(X) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())
    collapsed = True
    # TODO: Always use the svg_2d for the preview?
    preview = svg_2d((tuple((dim,) for dim in X.shape)))
    data_repr = short_data_repr_html(X)
    data_icon = _icon("icon-database")

    return (
        "<div class='ad-array-wrap'>"
        f"<input id='{data_id}' class='ad-array-in' type='checkbox' {collapsed}>"
        f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
        f"<div class='ad-array-preview ad-preview'><span>{preview}</span></div>"
        f"<div class='ad-array-data'>{data_repr}</div>"
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
        array_section(adata.X),
        collapsible_section(
            "layers",
            details=summarize_X(adata),
            n_items=len(adata.layers),
            enabled=True,
            collapsed=True,
        ),
        collapsible_section(
            "obs",
            details=summarize_obs(adata.obs),
            n_items=len(adata.obs.columns),
            enabled=True,
            collapsed=True,
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
