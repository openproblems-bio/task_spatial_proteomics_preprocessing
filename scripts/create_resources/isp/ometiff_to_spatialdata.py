#!/usr/bin/env python3
"""Convert a registered multiplexed OME-TIFF into a SpatialData Zarr store.

The ISP benchmark workflow takes a SpatialData Zarr store as its only input,
but nothing in the repository produces one: the published run used a local
conversion of the MCMICRO exemplar-001 dataset that was never committed. This
script closes that gap for any MCMICRO-style dataset.

Typical use, starting from the exemplar-001 download:

    # 1. register the raw cycles into one mosaic (what mcmicro's ashlar step does)
    ashlar raw/*.ome.tiff --ffp illumination/*-ffp.tif --dfp illumination/*-dfp.tif \
        --align-channel 0 -o exemplar-001-registered.ome.tif

    # 2. convert to the SpatialData store the workflow expects
    python ometiff_to_spatialdata.py \
        --image exemplar-001-registered.ome.tif \
        --markers exemplar-001/markers.csv \
        --output exemplar-001.zarr \
        --image-key exemplar-001_image

Optionally attach a segmentation mask as a labels element, so that the
`reference` segmentation method has something to read:

    ... --labels segmentation/cell.tif --labels-key exemplar-001_cell
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel


def read_markers(path: Path, n_channels: int) -> list[str]:
    """Read marker names from an MCMICRO-style markers.csv."""
    df = pd.read_csv(path)
    for col in ["marker_name", "Marker", "marker", "name"]:
        if col in df.columns:
            names = [str(x) for x in df[col].tolist()]
            break
    else:
        raise ValueError(
            f"No marker-name column in {path}. Columns present: {list(df.columns)}"
        )

    if len(names) != n_channels:
        raise ValueError(
            f"{path} lists {len(names)} markers but the image has {n_channels} channels."
        )
    return names


def read_image_cyx(path: Path) -> np.ndarray:
    """Read an OME-TIFF as a (channels, y, x) array."""
    with tifffile.TiffFile(path) as tif:
        arr = tif.series[0].asarray()

    if arr.ndim == 2:
        arr = arr[None, :, :]
    if arr.ndim != 3:
        raise ValueError(f"Expected a 2D or 3D image, got shape {arr.shape}.")

    # Registered mosaics are written channel-first, but be explicit about it
    # rather than trusting axis order silently: the channel axis is the short one.
    if arr.shape[0] > arr.shape[-1] and arr.shape[-1] <= 256:
        arr = np.moveaxis(arr, -1, 0)
    return arr


def read_labels_yx(path: Path) -> np.ndarray:
    with tifffile.TiffFile(path) as tif:
        arr = tif.series[0].asarray()
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D label mask, got shape {arr.shape}.")
    if arr.min() < 0:
        raise ValueError("Label masks must not contain negative values.")
    return arr.astype(np.int32, copy=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", required=True, type=Path, help="Registered multichannel OME-TIFF.")
    ap.add_argument("--markers", required=True, type=Path, help="markers.csv listing one marker per channel.")
    ap.add_argument("--output", required=True, type=Path, help="Destination .zarr store.")
    ap.add_argument("--image-key", default="image", help="Key for the image element.")
    ap.add_argument("--labels", type=Path, default=None, help="Optional segmentation mask TIFF.")
    ap.add_argument("--labels-key", default="cell", help="Key for the labels element.")
    args = ap.parse_args()

    print(f">> Reading {args.image}", flush=True)
    image = read_image_cyx(args.image)
    print(f"   shape {image.shape}, dtype {image.dtype}", flush=True)

    markers = read_markers(args.markers, image.shape[0])
    print(f">> {len(markers)} markers: {', '.join(markers)}", flush=True)

    elements = {
        "images": {
            args.image_key: Image2DModel.parse(image, dims=("c", "y", "x"), c_coords=markers)
        }
    }

    if args.labels is not None:
        print(f">> Reading labels from {args.labels}", flush=True)
        labels = read_labels_yx(args.labels)
        if labels.shape != image.shape[1:]:
            raise ValueError(
                f"Label shape {labels.shape} does not match image spatial shape {image.shape[1:]}."
            )
        n_objects = int((np.unique(labels) > 0).sum())
        print(f"   {n_objects} objects", flush=True)
        elements["labels"] = {args.labels_key: Labels2DModel.parse(labels, dims=("y", "x"))}

    print(f">> Writing {args.output}", flush=True)
    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; remove it first.")
    sdata = sd.SpatialData(**elements)
    sdata.write(args.output)

    print(">> Done", flush=True)
    print(sdata, flush=True)


if __name__ == "__main__":
    main()
