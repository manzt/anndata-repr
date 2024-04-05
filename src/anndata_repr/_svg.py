from __future__ import annotations

import typing

if typing.TYPE_CHECKING:
    import anndata

template = """
<svg xmlns="http://www.w3.org/2000/svg" height="300">

    <text x="140" y="15" font-size="1.0rem" font-weight="400" text-anchor="middle">11505</text>

    <g transform="translate(264, 42)">
        <text
            x="-35"
            y="10"
            dy="6"
            font-size="1rem"
            font-weight="400"
            text-anchor="middle"
            transform="rotate(-90, 0, 0)"
        >
            2638
        </text>
    </g>

    <g transform="translate(42, 20)">
        <!-- var #2B96C0 -->
        <rect x="0" y="0" width="200" rx="3" height="20" class="fill-sky-400 hover:fill-sky-500" />
        <text x="100" y="10" dy="4" font-size="0.8rem" font-weight="200" text-anchor="middle" fill="white">var</text>
    </g>

    <g transform="translate(42, 42)">
        <!-- X #50BA6F -->
        <rect x="0" y="0" width="200" height="70" class="fill-green-400 hover:fill-green-500" />
        <text x="100" y="35" dy="4" font-size="0.8rem" font-weight="200" text-anchor="middle" fill="white">X</text>
    </g>

    <g transform="translate(243, 42)">
        <!-- obs #EFC41B -->
        <rect x="0" y="0" rx="3" width="20" height="70" class="fill-yellow-400 hover:fill-yellow-500" />
        <text
            x="0"
            y="54"
            dy="-6"
            font-size="0.8rem"
            font-weight="200"
            text-anchor="middle"
            fill="white"
            transform="rotate(-90, 0, 35)"
        >
            obs
        </text>
    </g>

    <g transform="translate(21, 42)">
        <!-- obsm #EF9120 -->
        <rect x="0" y="0" rx="3" width="20" height="70" class="fill-orange-400 hover:fill-orange-500" />
        <text
            x="0"
            y="54"
            dy="-6"
            font-size="0.8rem"
            font-weight="200"
            text-anchor="middle"
            fill="white"
            transform="rotate(-90, 0, 35)"
        >
            obsm
        </text>
    </g>

    <g transform="translate(0, 42)">
        <!-- obsp #F15C5A -->
        <rect x="0" y="0" rx="3" width="20" height="70" class="fill-red-400 hover:fill-red-500" />
        <text
            x="0"
            y="54"
            dy="-6"
            font-size="0.8rem"
            font-weight="200"
            text-anchor="middle"
            fill="white"
            transform="rotate(-90, 0, 35)"
        >
            obsp
        </text>
    </g>

    <g transform="translate(42, 113)">
        <!-- varm #1A4C61 -->
        <rect x="0" y="0" rx="3" width="200" height="20" class="fill-sky-700 hover:fill-sky-900" />
        <text x="100" y="10" dy="4" font-size="0.8rem" font-weight="200" text-anchor="middle" fill="white">varm</text>
    </g>

    <g transform="translate(42, 134)">
        <!-- varp #965BA5 -->
        <rect x="0" y="0" rx="3" width="200" height="20" class="fill-purple-400 hover:fill-purple-500" />
        <text x="100" y="10" dy="4" font-size="0.8rem" font-weight="200" text-anchor="middle" fill="white">varp</text>
    </g>

</svg>
"""


def anndata_svg(adata: anndata.AnnData) -> str:
    shape = (adata.n_obs, adata.n_vars)
    return template
