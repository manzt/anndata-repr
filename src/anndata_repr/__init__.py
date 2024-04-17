from ._formatting_html import format_anndata_html  # noqa: F401
from ._version import __version__  # noqa: F401

try:
    import anndata
    anndata.AnnData._repr_html_ = format_anndata_html

except ImportError:
    import warnings

    warnings.warn("anndata not found. Please install anndata to enable anndata formatting.")
