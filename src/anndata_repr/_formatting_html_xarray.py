# This file has been modified from its original version in the xarray project,
# which is licensed under the Apache License, Version 2.0 (the "License").
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# The original file is located at: https://raw.githubusercontent.com/pydata/xarray/97d3a3aaa071fa5341132331abe90ec39f914b52/xarray/core/formatting_html.py
from __future__ import annotations

import contextlib
import math
import typing
import uuid
from datetime import datetime, timedelta
from functools import lru_cache
from html import escape
from importlib.resources import files
from itertools import chain, zip_longest
from pathlib import Path

import numpy as np
import pandas as pd


if typing.TYPE_CHECKING:
    duckarray = typing.Any

__all__ = [
    "_icon",
    "collapsible_section",
    "short_data_repr_html",
    "inline_variable_array_repr",
]


def last_n_items(array: pd.Series, n_desired: int):
    """Returns the last n_desired items of an array"""
    return np.ravel(np.asarray(array))[-n_desired:]


def first_n_items(array: pd.Series, n_desired: int):
    """Returns the first n_desired items of an array"""
    return np.ravel(np.asarray(array))[:n_desired]


def format_timedelta(t, timedelta_format=None):
    """Cast given object to a Timestamp and return a nicely formatted string"""
    try:
        timedelta_str = str(pd.Timedelta(t))
        days_str, time_str = timedelta_str.split(" days ")
    except Exception:
        # catch NaT and others that don't split nicely
        return str(t)
    else:
        if timedelta_format == "date":
            return days_str + " days"
        elif timedelta_format == "time":
            return time_str
        else:
            return timedelta_str


def format_timestamp(t):
    """Cast given object to a Timestamp and return a nicely formatted string"""
    try:
        timestamp = pd.Timestamp(t)
        datetime_str = timestamp.isoformat(sep=" ")
    except Exception:
        datetime_str = str(t)

    try:
        date_str, time_str = datetime_str.split()
    except ValueError:
        # catch NaT and others that don't split nicely
        return datetime_str
    else:
        if time_str == "00:00:00":
            return date_str
        else:
            return f"{date_str}T{time_str}"


def is_duck_array(
    value: typing.Any,
) -> typing.TypeGuard[duckarray[typing.Any, typing.Any]]:
    # TODO: replace is_duck_array with runtime checks via _arrayfunction_or_api protocol on
    # python 3.12 and higher (see https://github.com/pydata/xarray/issues/8696#issuecomment-1924588981)
    if isinstance(value, np.ndarray):
        return True
    return (
        hasattr(value, "ndim")
        and hasattr(value, "shape")
        and hasattr(value, "dtype")
        and (
            (hasattr(value, "__array_function__") and hasattr(value, "__array_ufunc__"))
            or hasattr(value, "__array_namespace__")
        )
    )


def short_array_repr(array: np.ndarray):
    # default to lower precision so a full (abbreviated) line can fit on
    # one line with the default display_width
    options = {
        "precision": 6,
        # "linewidth": OPTIONS["display_width"],
        # "threshold": OPTIONS["display_values_threshold"],
    }
    if array.ndim < 3:
        edgeitems = 3
    elif array.ndim == 3:
        edgeitems = 2
    else:
        edgeitems = 1
    options["edgeitems"] = edgeitems
    with set_numpy_options(**options):
        return repr(array)


@contextlib.contextmanager
def set_numpy_options(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally:
        np.set_printoptions(**original)


def limit_lines(string: str, *, limit: int):
    """
    If the string is more lines than the limit,
    this returns the middle lines replaced by an ellipsis
    """
    lines = string.splitlines()
    if len(lines) > limit:
        string = "\n".join(chain(lines[: limit // 2], ["..."], lines[-limit // 2 :]))
    return string


STATIC_FILES = (
    ("xarray.static.html", "icons-svg-inline.html"),
    ("xarray.static.css", "style.css"),
)


# @lru_cache(None)
def _load_static_files() -> tuple[str, str, str]:
    """Lazily load the resource files into memory the first time they are needed"""
    static_files = Path(__file__).parent / "static"
    return (
        (static_files / "icons-svg-inline.html").read_text(encoding="utf-8"),
        (static_files / "style.css").read_text(encoding="utf-8"),
        (static_files / "script.js").read_text(encoding="utf-8"),
    )


def short_data_repr_html(obj: pd.Series) -> str:
    if hasattr(obj, "_repr_html_"):
        return obj._repr_html_()
    text = short_array_repr(np.asarray(obj))
    return f"<pre>{escape(text)}</pre>"


def _icon(icon_name) -> str:
    # icon_name should be defined in xarray/static/html/icon-svg-inline.html
    return (
        f"<svg class='icon ad-{icon_name}'>"
        f"<use xlink:href='#{icon_name}'>"
        "</use>"
        "</svg>"
    )


def collapsible_section(
    name,
    inline_details: str = "",
    details: str = "",
    n_items: int | None = None,
    enabled: bool = True,
    collapsed: bool = False,
    color: str = "inherit",
) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-name" + str(uuid.uuid4())

    has_items = n_items is not None and n_items
    n_items_span = "" if n_items is None else f" <span>({n_items})</span>"
    enabled_str = "" if enabled and has_items else "disabled"
    collapsed_str = "" if collapsed or not has_items else "checked"
    tip = " title='Expand/collapse section'" if enabled else ""

    return (
        f"<input data-anndata='{name}' id='{data_id}' class='ad-section-summary-in'"
        f"type='checkbox' {enabled_str} {collapsed_str}>"
        f"<label for='{data_id}' class='ad-section-summary' {tip}>"
        f"{name}:{n_items_span}</label>"
        f"<div class='ad-section-inline-details'>{inline_details}</div>"
        f"<div class='ad-section-details'>{details}</div>"
    )


def _obj_repr(obj, header_components, sections):
    """Return HTML repr of an xarray object.

    If CSS is not injected (untrusted notebook), fallback to the plain text repr.

    """
    header = f"<div class='ad-header'>{''.join(h for h in header_components)}</div>"
    sections = "".join(f"<li class='ad-section-item'>{s}</li>" for s in sections)

    root_id = str(uuid.uuid4().hex)

    icons_svg, css_style, js_code = _load_static_files()
    js_code = js_code.replace("__ID__", root_id)
    return (
        f"<div id={root_id}>"
        f"{icons_svg}<style>{css_style}</style>"
        f'<script type="module">{js_code}</script>'
        f"<pre class='ad-text-repr-fallback'>{escape(repr(obj))}</pre>"
        "<div class='ad-wrap' style='display:none'>"
        f"{header}"
        f"<ul class='ad-sections'>{sections}</ul>"
        "</div>"
        "</div>"
    )


def format_item(x, timedelta_format=None, quote_strings=True):
    """Returns a succinct summary of an object as a string"""
    if isinstance(x, (np.datetime64, datetime)):
        return format_timestamp(x)
    if isinstance(x, (np.timedelta64, timedelta)):
        return format_timedelta(x, timedelta_format=timedelta_format)
    elif isinstance(x, (str, bytes)):
        if hasattr(x, "dtype"):
            x = x.item()  # type: ignore
        return repr(x) if quote_strings else x
    elif hasattr(x, "dtype") and np.issubdtype(x.dtype, np.floating):
        return f"{x.item():.4}"
    else:
        return str(x)


def format_items(x):
    """Returns a succinct summaries of all items in a sequence as strings"""
    x = np.asarray(x)
    timedelta_format = "datetime"
    if np.issubdtype(x.dtype, np.timedelta64):
        x = x.astype(dtype="timedelta64[ns]")
        day_part = x[~pd.isnull(x)].astype("timedelta64[D]").astype("timedelta64[ns]")
        time_needed = x[~pd.isnull(x)] != day_part
        day_needed = day_part != np.timedelta64(0, "ns")
        if np.logical_not(day_needed).all():
            timedelta_format = "time"
        elif np.logical_not(time_needed).all():
            timedelta_format = "date"

    formatted = [format_item(xi, timedelta_format) for xi in x]
    return formatted


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


def maybe_truncate(obj: typing.Any, maxlen: int = 500):
    s = str(obj)
    if len(s) > maxlen:
        s = s[: (maxlen - 3)] + "..."
    return s


def inline_variable_array_repr(col: pd.Series, max_width: int):
    """Build a one-line summary of a variable's data."""
    return format_array_flat(col, max_width)
