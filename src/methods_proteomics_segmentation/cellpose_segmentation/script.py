import numpy as np


def run(par: dict):
    import torch
    from cellpose import models

    model_type = par.get("model_type") or "cyto"

    image, nuc, mem, channel_info = load_image_and_channels(par, use_membrane=True)
    del image

    if mem is not None and model_type == "nuclei":
        # In two-channel mode Cellpose feeds the *first* channel to the model,
        # so the nuclei model would segment the membrane stain. Segment the
        # nuclear channel instead and record that the membrane was not used.
        print("WARNING: --model_type nuclei ignores --membrane_markers; segmenting the nuclear channel.", flush=True)
        mem = None
        channel_info["membrane"] = None

    if mem is not None:
        # Cellpose two-channel convention: channel 1 = cytoplasm/membrane,
        # channel 2 = nucleus (1-based), here stacked on axis 0.
        x = np.stack([mem, nuc]).astype(np.float32)
        eval_layout = {"channels": [1, 2], "channel_axis": 0}
    else:
        x = nuc.astype(np.float32)
        eval_layout = {"channels": [0, 0]}

    use_gpu = torch.cuda.is_available()
    print(f"Using {'GPU' if use_gpu else 'CPU'}", flush=True)

    model = models.Cellpose(gpu=use_gpu, model_type=model_type)

    eval_params = {
        k: par[k]
        for k in (
            "diameter", "batch_size", "resample", "normalize", "invert", "augment",
            "tile_overlap", "flow_threshold", "cellprob_threshold", "min_size", "niter",
        )
        if par.get(k) is not None
    }
    # Cellpose estimates the diameter with its size model when given 0/None.
    if eval_params.get("diameter") == 0:
        eval_params["diameter"] = None
    print(f"Running Cellpose ({model_type}) with {eval_layout} and {eval_params}", flush=True)
    masks, _, _, diams = model.eval(x, **eval_layout, **eval_params)

    labels = finalize_labels(masks)
    params = {"model_type": model_type, **eval_params, **eval_layout, "estimated_diameter": float(np.mean(diams))}
    write_outputs(par, labels, "cellpose_segmentation", params, channel_info)


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_markers": None,
    "nuclear_markers": None,
    "membrane_markers": None,
    "nuclei_channel": 0,
    "model_type": "cyto",
    "diameter": 30.0,
    "batch_size": 8,
    "resample": True,
    "normalize": True,
    "invert": False,
    "augment": False,
    "tile_overlap": 0.1,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0.0,
    "min_size": 15,
    "niter": 0,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "cellpose_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import finalize_labels, load_image_and_channels, write_outputs  # noqa: E402

if __name__ == "__main__":
    run(par)
