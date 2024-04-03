from ._version import __version__  # noqa: F401
from ._formatting_html import format_anndata_html  # noqa: F401


# just a hack to inject our _repr_html_ method into anndata objects
# TODO: remove with an offical integration with anndata
def read_h5ad(path: str):
    import anndata

    adata = anndata.read_h5ad(path)
    # just monkey patch the _repr_html_ method
    adata._repr_html_ = lambda: format_anndata_html(adata)  # type: ignore
    return adata
