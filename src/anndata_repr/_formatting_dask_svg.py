import math

import numpy as np

__all__ = ["svg_2d"]

text_style = 'font-size="1.0rem" font-weight="100" text-anchor="middle"'

def svg_lines(x1, y1, x2, y2, max_n=20):
    """Convert points into lines of text for an SVG plot

    Examples
    --------
    >>> svg_lines([0, 1], [0, 0], [10, 11], [1, 1])  # doctest: +NORMALIZE_WHITESPACE
    ['  <line x1="0" y1="0" x2="10" y2="1" style="stroke-width:2" />',
     '  <line x1="1" y1="0" x2="11" y2="1" style="stroke-width:2" />']
    """
    n = len(x1)

    if n > max_n:
        indices = np.linspace(0, n - 1, max_n, dtype="int")
    else:
        indices = range(n)

    lines = [
        '  <line x1="%d" y1="%d" x2="%d" y2="%d" />' % (x1[i], y1[i], x2[i], y2[i])
        for i in indices
    ]

    lines[0] = lines[0].replace(" /", ' style="stroke-width:2" /')
    lines[-1] = lines[-1].replace(" /", ' style="stroke-width:2" /')
    return lines


def svg_grid(x, y, offset=(0, 0), skew=(0, 0), size=200):
    """Create lines of SVG text that show a grid

    Parameters
    ----------
    x: numpy.ndarray
    y: numpy.ndarray
    offset: tuple
        translational displacement of the grid in SVG coordinates
    skew: tuple
    """
    # Horizontal lines
    x1 = np.zeros_like(y) + offset[0]
    y1 = y + offset[1]
    x2 = np.full_like(y, x[-1]) + offset[0]
    y2 = y + offset[1]

    if skew[0]:
        y2 += x.max() * skew[0]
    if skew[1]:
        x1 += skew[1] * y
        x2 += skew[1] * y

    min_x = min(x1.min(), x2.min())
    min_y = min(y1.min(), y2.min())
    max_x = max(x1.max(), x2.max())
    max_y = max(y1.max(), y2.max())
    max_n = size // 6

    h_lines = ["", "  <!-- Horizontal lines -->"] + svg_lines(x1, y1, x2, y2, max_n)

    # Vertical lines
    x1 = x + offset[0]
    y1 = np.zeros_like(x) + offset[1]
    x2 = x + offset[0]
    y2 = np.full_like(x, y[-1]) + offset[1]

    if skew[0]:
        y1 += skew[0] * x
        y2 += skew[0] * x
    if skew[1]:
        x2 += skew[1] * y.max()

    v_lines = ["", "  <!-- Vertical lines -->"] + svg_lines(x1, y1, x2, y2, max_n)

    # lightgreen
    color = "4fba6f" if len(x) < max_n and len(y) < max_n else "16a34a"
    # orange
    # color = "ECB172" if len(x) < max_n and len(y) < max_n else "8B4903"
    corners = f"{x1[0]},{y1[0]} {x1[-1]},{y1[-1]} {x2[-1]},{y2[-1]} {x2[0]},{y2[0]}"
    rect = [
        "",
        "  <!-- Colored Rectangle -->",
        f'  <polygon points="{corners}" style="fill:#{color}A0;stroke-width:0"/>',
    ]

    return h_lines + v_lines + rect, (min_x, max_x, min_y, max_y)


def draw_sizes(shape, size=200):
    """Get size in pixels for all dimensions"""
    mx = max(shape)
    ratios = [mx / max(0.1, d) for d in shape]
    ratios = [ratio_response(r) for r in ratios]
    return tuple(size / r for r in ratios)


def svg_2d(chunks, offset=(0, 0), skew=(0, 0), size=200, sizes=None):
    shape = tuple(map(sum, chunks))
    sizes = sizes or draw_sizes(shape, size=size)
    y, x = grid_points(chunks, sizes)

    lines, (min_x, max_x, min_y, max_y) = svg_grid(
        x, y, offset=offset, skew=skew, size=size
    )

    header = (
        '<svg width="%d" height="%d" style="stroke:rgb(0,0,0);stroke-width:1" >\n'
        % (max_x + 50, max_y + 50)
    )
    footer = "\n</svg>"

    if shape[0] >= 100:
        rotate = -90
    else:
        rotate = 0

    text = [
        "",
        "  <!-- Text -->",
        f'  <text x="{max_x / 2}" y="{max_y + 20}" {text_style} >{shape[1]:,}</text>',
        f'  <text x="{max_x + 20}" y="{max_y / 2}" {text_style} transform="rotate({rotate},{max_x + 20},{max_y / 2})">{shape[0]:,}</text>',
    ]

    return header + "\n".join(lines + text) + footer


def ratio_response(x):
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


def grid_points(chunks, sizes):
    cumchunks = [np.cumsum((0,) + c) for c in chunks]
    points = [x * size / x[-1] for x, size in zip(cumchunks, sizes)]
    return points


