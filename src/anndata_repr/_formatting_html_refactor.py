from __future__ import annotations

import typing
import uuid
from html import escape
import math
from datetime import datetime, timedelta

import numpy as np

from ._formatting_dask_svg import svg_2d
from ._formatting_html_xarray import (
    _icon,
    _obj_repr,
    collapsible_section,
    short_data_repr_html,
    inline_variable_array_repr,
    format_timestamp,
    format_timedelta
)

if typing.TYPE_CHECKING:
    import anndata
    import pandas as pd

__all__ = ["format_anndata_html"]

def summarize_attrs(obj: pd.Series) -> dict:
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
    return enum

def summarize_columns(name: str, col: pd.Series, is_index: bool = True) -> dict:
    """Summarize a single column of a DataFrame.

    Parameters
    ----------
    name : str
        The name of the column.

    col : pd.Series
        The column to summarize.

    Returns
    -------
    dict
        Dictionary with information
    """
    name = escape(str(name))
    dtype = escape(str(col.dtype))

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())

    preview = escape(inline_variable_array_repr(col, 35))
    attrs_ul = summarize_attrs(col)
    data_repr = short_data_repr_html(col)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    obj = {
        'name': name,
        'dtype': dtype,
        'preview': preview,
        'attrs': {
            'id': attrs_id,
            'icon': attrs_icon,
            'content': attrs_ul,
        },
        'data': {
            'id': data_id,
            'icon': data_icon,
            'content': data_repr,
        }
    }

    return obj



def summarize_table(df: pd.DataFrame, is_index: bool = True) -> dict:

    # "unique" ids required to expand/collapse subsections
    data_id = "data-" + str(uuid.uuid4())

    data_icon = _icon("icon-database")

    obj = {
        "id": data_id,
        "icon": data_icon,
        "content": df[0:10]
    }

    return obj


def summarize_df(df: pd.DataFrame) -> dict:
    columns_list = []
    for k in df: 
        assert isinstance(k, str), "Column of dataframe is not a string"
        columns_list.append(summarize_columns(k, df[k]))
    
    obj = {
        "table": summarize_table(df),
        "columns": columns_list
    }
    return obj


# def format_var_obs(adata: anndata.AnnData) -> str:
#     dims_li = "".join(
#         f"<li><span class='xr-has-index'>obs</span>: {adata.n_obs}</li>"
#         f"<li><span class='xr-has-index'>var</span>: {adata.n_vars}</li>"
#     )
#     return f"<ul class='xr-dim-list'>{dims_li}</ul>"


# def array_section(X) -> str:
#     # "unique" id to expand/collapse the section
#     data_id = "section-" + str(uuid.uuid4())
#     collapsed = True
#     # TODO: Always use the svg_2d for the preview?
#     preview = svg_2d((tuple((dim,) for dim in X.shape)))
#     data_repr = short_data_repr_html(X)
#     data_icon = _icon("icon-database")

#     return (
#         "<div class='xr-array-wrap'>"
#         f"<input id='{data_id}' class='xr-array-in' type='checkbox' {collapsed}>"
#         f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
#         f"<div class='xr-array-preview xr-preview'><span>{preview}</span></div>"
#         f"<div class='xr-array-data'>{data_repr}</div>"
#         "</div>"
#     )


def format_anndata_html(adata: anndata.AnnData) -> str:
    """Format an AnnData object as an HTML string.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to format.

    """
    obj_type = f"anndata.{type(adata).__name__}"
    arr_name = ""  # TODO: add somethign here?

    # header_components = [
    #     f"<div class='xr-obj-type'>{obj_type}</div>",
    #     f"<div class='xr-array-name'>{arr_name}</div>",
    #     format_var_obs(adata),
    # ]

    return str(summarize_df(adata.obs))
    # print(summarize_df(adata.obs))

    # sections = [
    #     array_section(adata.X),
    #     collapsible_section(
    #         "obs",
    #         details=summarize_df(adata.obs),
    #         n_items=len(adata.obs.columns),
    #         enabled=True,
    #         collapsed=False,
    #     ),
    #     collapsible_section(
    #         "var",
    #         details=summarize_df(adata.var),
    #         n_items=len(adata.var.columns),
    #         enabled=True,
    #         collapsed=True,
    #     ),
    # ]

    # return _obj_repr(adata, header_components, sections)





# Series


# def last_n_items(array: pd.Series, n_desired: int):
#     """Returns the last n_desired items of an array"""
#     return np.ravel(np.asarray(array))[-n_desired:]


# def first_n_items(array: pd.Series, n_desired: int):
#     """Returns the first n_desired items of an array"""
#     return np.ravel(np.asarray(array))[:n_desired]


# def format_item(x, timedelta_format=None, quote_strings=True):
#     """Returns a succinct summary of an object as a string"""
#     if isinstance(x, (np.datetime64, datetime)):
#         return format_timestamp(x)
#     if isinstance(x, (np.timedelta64, timedelta)):
#         return format_timedelta(x, timedelta_format=timedelta_format)
#     elif isinstance(x, (str, bytes)):
#         if hasattr(x, "dtype"):
#             x = x.item()  # type: ignore
#         return repr(x) if quote_strings else x
#     elif hasattr(x, "dtype") and np.issubdtype(x.dtype, np.floating):
#         return f"{x.item():.4}"
#     else:
#         return str(x)


# def format_items(x):
#     """Returns a succinct summaries of all items in a sequence as strings"""
#     x = np.asarray(x)
#     timedelta_format = "datetime"
#     if np.issubdtype(x.dtype, np.timedelta64):
#         x = x.astype(dtype="timedelta64[ns]")
#         day_part = x[~pd.isnull(x)].astype("timedelta64[D]").astype("timedelta64[ns]")
#         time_needed = x[~pd.isnull(x)] != day_part
#         day_needed = day_part != np.timedelta64(0, "ns")
#         if np.logical_not(day_needed).all():
#             timedelta_format = "time"
#         elif np.logical_not(time_needed).all():
#             timedelta_format = "date"

#     formatted = [format_item(xi, timedelta_format) for xi in x]
#     return formatted


# # from xarray.core.formatting import format_array_flat
# def format_array_flat(array: pd.Series, max_width: int):
#     """Return a formatted string for as many items in the flattened version of
#     array that will fit within max_width characters.
#     """
#     # every item will take up at least two characters, but we always want to
#     # print at least first and last items
#     max_possibly_relevant = min(max(array.size, 1), max(math.ceil(max_width / 2.0), 2))
#     relevant_front_items = format_items(
#         first_n_items(array, (max_possibly_relevant + 1) // 2)
#     )
#     relevant_back_items = format_items(last_n_items(array, max_possibly_relevant // 2))
#     # interleave relevant front and back items:
#     #     [a, b, c] and [y, z] -> [a, z, b, y, c]
#     relevant_items = sum(
#         zip_longest(relevant_front_items, reversed(relevant_back_items)), ()
#     )[:max_possibly_relevant]

#     cum_len = np.cumsum([len(s) + 1 for s in relevant_items]) - 1
#     if (array.size > 2) and (
#         (max_possibly_relevant < array.size) or (cum_len > max_width).any()
#     ):
#         padding = " ... "
#         max_len = max(int(np.argmax(cum_len + len(padding) - 1 > max_width)), 2)
#         count = min(array.size, max_len)
#     else:
#         count = array.size
#         padding = "" if (count <= 1) else " "

#     num_front = (count + 1) // 2
#     num_back = count - num_front
#     # note that num_back is 0 <--> array.size is 0 or 1
#     #                         <--> relevant_back_items is []
#     pprint_str = "".join(
#         [
#             " ".join(relevant_front_items[:num_front]),  # type: ignore
#             padding,
#             " ".join(relevant_back_items[-num_back:]),  # type: ignore
#         ]
#     )

#     # As a final check, if it's still too long even with the limit in values,
#     # replace the end with an ellipsis
#     # NB: this will still returns a full 3-character ellipsis when max_width < 3
#     if len(pprint_str) > max_width:
#         pprint_str = pprint_str[: max(max_width - 3, 0)] + "..."

#     return pprint_str




# render
# def dataframe_to_table(dataframe):
#     # Initialize the table with the correct class names for styling
#     table = """<div class='relative overflow-x-auto my-div'>
#     <table class="">
#     """

#     def make_header(columns):
#         # Style the header row according to the provided CSS class names
#         header = '<thead class="table-header"><tr>'
#         for column in columns:
#             header += '<th scope="col" class="column-header">'+column+"</th>"
#         header += "</tr></thead>"
#         return header

#     def make_truncated_data(data):
#         # Function to truncate data if it's longer than 10 characters
#         if len(str(data)) > 10:
#             return str(data)[:10] + "..."
#         return str(data)

#     def make_row(row, index):
#         # Style each row according to the provided CSS class names, alternating row color not implemented in CSS
#         row_html = f'<tr class="table-row">'
#         for value in row:
#             row_html += f'<td class="table-cell">{make_truncated_data(value)}</td>'
#         row_html += "</tr>"
#         return row_html

#     # Construct the table header
#     table += make_header(dataframe.columns)

#     # Construct each row of the table
#     for index, row in enumerate(dataframe.itertuples(index=False), start=1):
#         table += make_row(row, index)

#     # Close the table and div tags
#     table += "</table></div>"
#     return table