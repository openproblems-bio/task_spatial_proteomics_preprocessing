"""Tiny synthetic fixtures for ISP component smoke tests.

The ISP benchmark was developed against a single large MCMICRO dataset that is
not redistributable and takes minutes to process. These helpers build
deliberately small, deterministic stand-ins (64x64 pixels, 3 channels) so that
every component can be exercised in CI without any external resource.

Imports are kept local to each helper: component images only install the
dependencies that component needs, so a segmentation test must not pay for
anndata and a QC test must not pay for spatialdata.
"""

IMAGE_SHAPE = (4, 64, 64)
# Two nuclear replicates (the same nuclei imaged in two cycles, as multiplexed
# panels do) plus a mutually exclusive immune/epithelial pair. This shape is
# what lets the segmentation-quality metrics be tested at all.
MARKERS = ["DNA_1", "DNA_2", "CD45", "ECAD"]
NUCLEAR_MARKERS = ["DNA_1", "DNA_2"]
EXCLUSIVE_PAIR = ("CD45", "ECAD")

# Blob centres used by both the image and the label fixtures, with a deliberate
# range of radii so that min_area filtering has something to remove.
_BLOBS = [(14, 14, 5), (14, 46, 5), (46, 14, 5), (46, 46, 4), (30, 30, 2)]


def build_image(seed: int = 0):
    """Return a (3, 64, 64) float32 channel-y-x image with bright nuclear blobs."""
    import numpy as np

    rng = np.random.default_rng(seed)
    img = (rng.random(IMAGE_SHAPE) * 10.0).astype(np.float32)

    yy, xx = np.ogrid[: IMAGE_SHAPE[1], : IMAGE_SHAPE[2]]
    for cy, cx, radius in _BLOBS:
        disk = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
        # Both nuclear channels light up the same blobs, as replicate stains of
        # the same nuclei do, with a little independent noise so their
        # correlation is high but not exactly 1.
        img[0][disk] += 500.0
        img[1][disk] += 500.0 + rng.normal(0.0, 5.0, size=disk.sum()).astype(np.float32)
        img[2][disk] += 120.0  # CD45, tracks the nuclei
    # The last channel gets a smooth left-to-right illumination gradient so that
    # flat-field correction has a low-frequency field to remove.
    img[3] += np.linspace(0.0, 200.0, IMAGE_SHAPE[2], dtype=np.float32)[None, :]

    return img


def build_labels():
    """Return a (64, 64) int32 label mask with gapless labels 1..n."""
    import numpy as np

    labels = np.zeros(IMAGE_SHAPE[1:], dtype=np.int32)
    yy, xx = np.ogrid[: IMAGE_SHAPE[1], : IMAGE_SHAPE[2]]
    for i, (cy, cx, radius) in enumerate(_BLOBS, start=1):
        labels[(yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2] = i
    return labels


def make_image(path, seed: int = 0):
    """Write the fixture image to a .npy file."""
    import numpy as np

    img = build_image(seed=seed)
    np.save(path, img)
    return img


def make_labels(path):
    """Write the fixture label mask to a .npy file."""
    import numpy as np

    labels = build_labels()
    np.save(path, labels)
    return labels


def make_markers_json(path, markers=None):
    """Write the marker sidecar produced by the image-correction components."""
    import json

    markers = list(markers) if markers is not None else list(MARKERS)
    with open(path, "w") as f:
        json.dump({"method": "fixture", "image_key": "fixture_image", "markers": markers}, f)
    return markers


def make_cell_by_protein(
    path,
    n_cells: int = 40,
    seed: int = 0,
    with_normalized: bool = False,
    with_qc: bool = False,
    with_celltypes: bool = False,
):
    """Write a cell-by-protein AnnData matching the aggregation output contract."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    X = (rng.lognormal(mean=3.0, sigma=1.0, size=(n_cells, len(MARKERS)))).astype(np.float32)

    # DNA_2 is the same nuclei restained, so it tracks DNA_1 closely. Leaving
    # them independent would make replicate concordance read ~0 on data that is
    # by construction fine.
    i0, i1 = MARKERS.index("DNA_1"), MARKERS.index("DNA_2")
    X[:, i1] = X[:, i0] * rng.normal(1.0, 0.05, size=n_cells).astype(np.float32)
    X[:, i1] = np.abs(X[:, i1])

    # Applied before anything is derived from X, so counts, the normalized
    # layer and the QC totals all describe the same matrix.
    cell_types = None
    if with_celltypes:
        # Two groups that differ in marker *composition*, not just in overall
        # brightness. A uniform intensity shift would be largely undone by
        # per-cell total normalization, which is the layer the separability
        # metric scores, so the groups would not actually be separable there.
        group = rng.integers(0, 2, size=n_cells)
        X[group == 1, MARKERS.index("CD45")] *= 20.0
        X[group == 0, MARKERS.index("ECAD")] *= 20.0
        cell_types = np.where(group == 1, "type_a", "type_b")

    # A quarter of the cells are debris: small *and* dim, which is what QC
    # exists to catch. Making this a fixed fraction rather than a lucky draw
    # keeps QC-sensitive assertions deterministic - with areas drawn uniformly
    # the number falling below a min_area of 50 varied with the RNG stream, and
    # at one point only 2 of 40 cells were rejected.
    n_debris = max(1, n_cells // 4)
    is_debris = np.zeros(n_cells, dtype=bool)
    is_debris[:n_debris] = True
    X[is_debris] *= 0.05

    area = np.empty(n_cells, dtype=float)
    area[is_debris] = rng.integers(5, 45, size=n_debris)
    area[~is_debris] = rng.integers(60, 400, size=n_cells - n_debris)

    obs = pd.DataFrame(
        {
            "cell_id": [str(i + 1) for i in range(n_cells)],
            "area": area,
            "is_debris": is_debris,
            "centroid_y": rng.random(n_cells) * 64.0,
            "centroid_x": rng.random(n_cells) * 64.0,
        },
        index=[str(i + 1) for i in range(n_cells)],
    )

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=list(MARKERS)))
    adata.layers["counts"] = X.copy()
    if cell_types is not None:
        adata.obs["cell_type"] = cell_types
    if with_qc:
        # Stand in for the QC step: annotate, do not filter.
        adata.obs["passed_QC"] = (adata.obs["area"] >= 50).to_numpy()
        adata.obs["total_intensity"] = X.sum(axis=1).astype(float)
    if with_normalized:
        totals = X.sum(axis=1, keepdims=True)
        totals[totals <= 0] = 1.0
        adata.layers["normalized"] = np.log1p(X / totals * 1e4).astype(np.float32)
    adata.uns["aggregation"] = {"method": "fixture", "statistic": "mean"}

    adata.write_h5ad(path)
    return adata


def make_spatialdata_zarr(path, seed: int = 0):
    """Write a minimal SpatialData store with one image and one labels element."""
    import spatialdata as sd
    from spatialdata.models import Image2DModel, Labels2DModel

    img = build_image(seed=seed)
    labels = build_labels()

    sdata = sd.SpatialData(
        images={
            "fixture_image": Image2DModel.parse(
                img, dims=("c", "y", "x"), c_coords=list(MARKERS)
            )
        },
        labels={"fixture_cell": Labels2DModel.parse(labels, dims=("y", "x"))},
    )
    sdata.write(path)
    return sdata
