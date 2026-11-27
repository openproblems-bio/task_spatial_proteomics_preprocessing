import numpy as np
import anndata as ad

try:
    import scipy.sparse as sp
except Exception:  # pragma: no cover
    sp = None


def _dense(x) -> np.ndarray:
    if sp is not None and sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def _get_counts_matrix(adata: ad.AnnData) -> np.ndarray:
    if "counts" in adata.layers:
        return _dense(adata.layers["counts"])
    return _dense(adata.X)


def run(par: dict):
    in_path = par["input"]
    out_path = par["output"]
    cofactor = float(par.get("cofactor", 5.0))

    if cofactor <= 0:
        raise ValueError("--cofactor must be > 0")

    adata = ad.read_h5ad(in_path)
    X = _get_counts_matrix(adata).astype(np.float32)

    Xn = np.arcsinh(X / cofactor).astype(np.float32)

    adata.layers["normalized"] = Xn
    adata.uns.setdefault("normalization", {})
    adata.uns["normalization"].update(
        {
            "method": "arcsinh_normalization",
            "description": "Cytometry-style arcsinh(x / cofactor) transformation.",
            "scale": None,
            "cofactor": cofactor,
        }
    )

    adata.write_h5ad(out_path, compression="gzip")


## VIASH START
par = {
    "input": "input.h5ad",
    "scale": 10000.0,
    "cofactor": 5.0,
    "transform": "log1p",
    "output": "output.h5ad",
}
meta = {"name": "arcsinh_normalization"}
## VIASH END


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        import argparse

        ap = argparse.ArgumentParser()
        ap.add_argument("--input", required=True)
        ap.add_argument("--scale", type=float, default=10000.0)
        ap.add_argument("--cofactor", type=float, default=5.0)
        ap.add_argument("--transform", default="log1p")
        ap.add_argument("--output", required=True)
        run(vars(ap.parse_args()))
