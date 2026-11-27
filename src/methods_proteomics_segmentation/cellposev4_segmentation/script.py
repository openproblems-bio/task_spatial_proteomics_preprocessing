import numpy as np


def run(par: dict):
    import torch
    from cellpose.models import CellposeModel

    image, nuc, mem, channel_info = load_image_and_channels(par, use_membrane=True)
    del image

    if mem is not None:
        x = np.stack([mem, nuc]).astype(np.float32)
        layout = {"channel_axis": 0}
    else:
        x = nuc.astype(np.float32)
        layout = {}

    use_gpu = torch.cuda.is_available()
    print(f"Using {'GPU' if use_gpu else 'CPU'}", flush=True)
    model = CellposeModel(gpu=use_gpu)

    eval_params = {
        k: par[k]
        for k in ("diameter", "flow_threshold", "cellprob_threshold", "niter", "min_size", "resample")
        if par.get(k) is not None
    }
    if eval_params.get("diameter") == 0:
        eval_params["diameter"] = None
    print(f"Running Cellpose-SAM with {layout} and {eval_params}", flush=True)
    masks, _, _ = model.eval(x, progress=True, **layout, **eval_params)

    labels = finalize_labels(masks)
    write_outputs(par, labels, "cellposev4_segmentation", {"model": "cpsam", **eval_params}, channel_info)


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_markers": None,
    "nuclear_markers": None,
    "membrane_markers": None,
    "nuclei_channel": 0,
    "diameter": 30.0,
    "flow_threshold": 0.0,
    "cellprob_threshold": 0.0,
    "niter": 10,
    "min_size": -1,
    "resample": False,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "cellposev4_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import finalize_labels, load_image_and_channels, write_outputs  # noqa: E402

if __name__ == "__main__":
    run(par)
