import numpy as np

# numpy>=1.24 removed the deprecated scalar aliases (np.bool, np.int, ...) that
# stardist/csbdeep still reference; restore them as numpy itself defined them.
for _alias, _target in {
    "bool": bool, "int": int, "float": float, "complex": complex,
    "object": object, "str": str, "long": int, "unicode": str,
}.items():
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)

BIG_PX = 4096


def _n_tiles(px):
    # Forward-pass tiling to bound GPU memory (~2048 px/tile); no stitching constraint.
    n = max(1, int(px) // 2048)
    return (n, n)


def segment_stardist(nuc, par):
    from csbdeep.data import Normalizer, normalize_mi_ma
    from stardist.models import StarDist2D

    model = StarDist2D.from_pretrained(par.get("model") or "2D_versatile_fluo")

    # Normalize with global percentiles rather than per tile, so tiles of a big
    # image see the same intensity scale (stardist predict_big_data example).
    class GlobalNormalizer(Normalizer):
        def __init__(self, mi, ma):
            self.mi, self.ma = mi, ma

        def before(self, x, axes):
            return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)

        def after(*args, **kwargs):
            assert False

        @property
        def do_after(self):
            return False

    finite = nuc[np.isfinite(nuc)]
    mi, ma = np.percentile(finite, [1, 99.8])
    normalizer = GlobalNormalizer(mi, ma)
    nuc = np.nan_to_num(nuc, nan=float(mi))

    eval_params = {k: par[k] for k in ("prob_thresh", "nms_thresh", "scale") if par.get(k) is not None}
    print(f"stardist overrides: {eval_params}", flush=True)

    # Images up to BIG_PX take the single-pass path, which has no stitching
    # constraint on object size. Larger ones are tiled and stitched, which
    # requires every object to be smaller than min_overlap and block_size to
    # comfortably exceed min_overlap + 2*context.
    max_dim = int(max(nuc.shape))
    if max_dim <= BIG_PX:
        labels, _ = model.predict_instances(
            nuc, axes="YX", normalizer=normalizer, n_tiles=_n_tiles(max_dim), **eval_params
        )
        return labels, eval_params

    block_size = BIG_PX
    min_overlap = int(par.get("max_object_diameter") or 192)
    context = int(min(min_overlap, 128))
    while True:
        try:
            labels, _ = model.predict_instances_big(
                nuc, axes="YX", block_size=block_size, min_overlap=min_overlap, context=context,
                n_tiles=_n_tiles(block_size), normalizer=normalizer, **eval_params,
            )
            return labels, {**eval_params, "min_overlap": min_overlap, "block_size": block_size}
        except RuntimeError as e:
            # A rare oversized blob exceeded min_overlap: double it and retry.
            if "min_overlap" not in str(e) or min_overlap >= 2048:
                raise
            min_overlap *= 2
            block_size = max(block_size, min_overlap + 2 * context + 256)
            print(
                "predict_instances_big: an object exceeded min_overlap; retrying with "
                f"min_overlap={min_overlap}, block_size={block_size}",
                flush=True,
            )


def run(par: dict):
    image, nuc, _, channel_info = load_image_and_channels(par)
    del image

    labels, used = segment_stardist(nuc, par)
    labels = finalize_labels(labels)
    write_outputs(
        par, labels, "stardist_segmentation",
        {"model": par.get("model") or "2D_versatile_fluo", **used}, channel_info,
    )


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_markers": None,
    "nuclear_markers": None,
    "nuclei_channel": 0,
    "model": "2D_versatile_fluo",
    "prob_thresh": None,
    "nms_thresh": None,
    "scale": None,
    "max_object_diameter": None,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "stardist_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import finalize_labels, load_image_and_channels, write_outputs  # noqa: E402

if __name__ == "__main__":
    run(par)
