from __future__ import annotations

from dataclasses import dataclass
import typing
import math

if typing.TYPE_CHECKING:
    import anndata


@dataclass
class Color:
    primary: str

    @property
    def hover(self):
        return f"color-mix(in srgb, {self.primary}, white 20%)"


_COLORS = {
    "X": Color("#3c8b53"),
    "var": Color("#2c96c0"),
    "obs": Color("#efc41c"),
    "obsm": Color("#ef9021"),
    "obsp": Color("#f15c5a"),
    "varm": Color("#194c61"),
    "varp": Color("#965ba5"),
}


def ratio_response(x: float):
    """How we display actual size ratios

    Common ratios in sizes span several orders of magnitude,
    which is hard for us to perceive.

    We keep ratios in the 1-3 range accurate, and then apply a logarithm to
    values up until about 100 or so, at which point we stop scaling.
    """
    if x < math.e:
        return x
    elif x <= 100:
        return math.log(x + 12.4)  # f(e) == e
    else:
        return math.log(100 + 12.4)


def draw_sizes(shape: tuple[int, ...], size: int):
    """Get size in pixels for all dimensions"""
    mx = max(shape)
    ratios = [mx / max(0.1, d) for d in shape]
    ratios = [ratio_response(r) for r in ratios]
    return tuple(size / r for r in ratios)


def style_tag():
    def hover_colors():
        for name, color in _COLORS.items():
            yield f"""
            #{name}:hover {{
                fill: {color.hover};
            }}
            """

    styles = [
        "<style>",
        "text { pointer-events: none; }",
        "rect { cursor: pointer; }",
        *hover_colors(),
        "</style>",
    ]
    return "\n".join(styles)


def anndata_svg(
    adata: anndata.AnnData, size: int = 200, spacing: int = 1, handle_size: int = 15
) -> str:
    h, w = draw_sizes(shape=(adata.n_obs, adata.n_vars), size=size)
    total_width = w + (handle_size + spacing) * 4
    total_height = h + (handle_size + spacing) * 4

    dim_text_style = 'font-size="1.0em" font-weight="400" text-anchor="middle"'
    dim_text_style_h = f'{dim_text_style} dy="0.4em"'
    dim_text_style_v = f'{dim_text_style} dy="0.4em" transform="rotate(-90, {handle_size / 2 }, {h / 2})"'

    handle_text_style = (
        'font-size="0.8em" font-weight="200" text-anchor="middle" fill="white"'
    )
    handle_text_style_h = f'{handle_text_style} dy="0.4em"'
    handle_text_style_v = f'{handle_text_style} dy="0.3em" transform="rotate(-90, {handle_size / 2}, {h / 2})"'
    rx = 0  # corner radius

    svg = [
        f'<svg class="ad-svg" xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{total_height}" shape-rendering="crispEdges" text-rendering="geometricPrecision" tabindex="0">',
        "<!-- n_vars -->",
        f'<g transform="translate({(handle_size + spacing) * 2}, 0)">'
        f'<text x="{w / 2}" y="{handle_size / 2}" {dim_text_style_h}>{adata.n_vars}</text>',
        "</g>",
        "<!-- var -->",
        f'<g transform="translate({(handle_size + spacing) * 2}, {handle_size + spacing})" role="button">'
        f'<rect id="var" width="{w}" height="{handle_size}" rx={rx} fill="{_COLORS["var"].primary}" />',
        f'<text x="{w / 2}" y="{handle_size / 2}" {handle_text_style_h}>var</text>',
        "</g>",
        "<!-- obsp -->",
        f'<g transform="translate(0, {(handle_size + spacing) * 2})" role="button">'
        f'<rect id="obsp" width="{handle_size}" height="{h}" rx={rx} fill="{_COLORS["obsp"].primary}" />',
        f'<text x="{handle_size / 2}" y="{h / 2}" {handle_text_style_v}>obsp</text>',
        "</g>",
        "<!-- obsm -->",
        f'<g transform="translate({handle_size + spacing}, {(handle_size + spacing) * 2})" role="button">'
        f'<rect id="obsm" width="{handle_size}" height="{h}" rx={rx} fill="{_COLORS["obsm"].primary}" />',
        f'<text x="{handle_size / 2}" y="{h / 2}" {handle_text_style_v}>obsm</text>',
        "</g>",
        "<!-- X -->",
        f'<g transform="translate({(handle_size + spacing) * 2}, {(handle_size + spacing) * 2})" role="button">'
        f'<rect id="X" width="{w}" height="{h}" rx={rx} fill="{_COLORS["X"].primary}" />',
        f'<text x="{w / 2}" y="{h / 2}" {handle_text_style_h}>X</text>',
        "</g>",
        "<!-- obs -->",
        f'<g transform="translate({(handle_size + spacing) * 2 + w + spacing}, {(handle_size + spacing) * 2})" role="button">'
        f'<rect id="obs" width="{handle_size}" height="{h}" rx={rx} fill="{_COLORS["obs"].primary}" />',
        f'<text x="{handle_size / 2}" y="{h / 2}" {handle_text_style_v}>obs</text>',
        "</g>",
        "<!-- varm -->",
        f'<g transform="translate({(handle_size + spacing) * 2}, {(handle_size + spacing) * 2 + (h + spacing)})" role="button">'
        f'<rect id="varm" width="{w}" height="{handle_size}" rx={rx} fill="{_COLORS["varm"].primary}" />',
        f'<text x="{w / 2}" y="{handle_size / 2}" {handle_text_style_h}>varm</text>',
        "</g>",
        "<!-- varp -->",
        f'<g transform="translate({(handle_size + spacing) * 2}, {(handle_size + spacing) * 3 + (h + spacing)})" role="button">'
        f'<rect id="varp" width="{w}" height="{handle_size}" rx={rx} fill="{_COLORS["varp"].primary}" />',
        f'<text x="{w / 2}" y="{handle_size / 2}" {handle_text_style_h}>varp</text>',
        "</g>",
        "<!-- uns -->",
        f'<g transform="translate(0, {(handle_size + spacing) * 2 + (h + spacing)})" role="button">'
        f'<rect id="uns" width="{(handle_size * 2) + spacing}" height="{(handle_size * 2) + spacing}" rx={rx} fill="transparent" />',
        f'<text x="{handle_size + spacing}" y="{handle_size + spacing}" dx="-0.1em" dy="0.25em" {dim_text_style}>{{uns}}</text>'
        "</g>",
        "<!-- n_obs -->",
        f'<g transform="translate({(handle_size + spacing) * 3 + w + spacing}, {(handle_size + spacing) * 2})">'
        f'<text x="{handle_size / 2}" y="{h / 2}" {dim_text_style_v}>{adata.n_obs}</text>',
        "</g>",
        style_tag(),
        "</svg>",
    ]

    return "\n".join(svg)
