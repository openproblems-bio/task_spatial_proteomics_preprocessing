#!/usr/bin/env bash

# Turn the public MCMICRO exemplar-001 download into the SpatialData Zarr store
# that the ISP benchmark workflow expects.
#
# The workflow takes a .zarr as its only input, but the published exemplar-001
# download is *raw* acquisition data: unstitched tiles for three imaging cycles,
# plus illumination profiles. Two steps bridge the gap:
#
#   1. registration - stitch and align the cycles into one 12-channel mosaic,
#      with illumination correction. This is what mcmicro's registration step
#      does, using the same tool (ashlar).
#   2. conversion   - wrap the mosaic and its marker names in a SpatialData store.
#
# Usage:
#   bash scripts/create_resources/isp/prepare_exemplar_001.sh [OUTPUT_DIR]
#
# The resulting store is what you pass as ZARR_PATH:
#   ZARR_PATH=<OUTPUT_DIR>/exemplar-001.zarr \
#   bash src/workflows/run_benchmark/test.sh
#
# Requirements: curl, unzip, and a Python environment with ashlar, tifffile,
# spatialdata, numpy and pandas. ashlar reads OME-TIFF through BioFormats, so a
# real JDK must be on PATH (javac included, not just a JRE).
#
# Note on reference segmentation: the public download contains no segmentation
# masks, so the `reference` segmentation method has nothing to read. Those
# labels are mcmicro segmentation output. If you have them, attach them with
# --labels/--labels-key when running ometiff_to_spatialdata.py, and pass the key
# to the workflow as LABELS_KEY.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:-${OUTPUT_DIR:-resources/isp}}"
URL="https://mcmicro.s3.amazonaws.com/exemplars/exemplar-001.zip"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo ">> Output directory: $(pwd)"

if ! command -v javac >/dev/null 2>&1; then
  echo "ERROR: no javac on PATH. ashlar reads OME-TIFF via BioFormats and needs a" >&2
  echo "       full JDK, not just a JRE. Install one and set JAVA_HOME." >&2
  exit 1
fi

if [ ! -f exemplar-001.zip ]; then
  echo ">> Downloading exemplar-001 (~240 MB)"
  curl -fL --retry 3 --retry-delay 5 -o exemplar-001.zip.part "$URL"
  mv exemplar-001.zip.part exemplar-001.zip
else
  echo ">> exemplar-001.zip already present, skipping download"
fi

if [ ! -d exemplar-001 ]; then
  echo ">> Extracting"
  unzip -q exemplar-001.zip
fi

if [ ! -f exemplar-001-registered.ome.tif ]; then
  echo ">> Registering the three cycles into one mosaic (this takes a few minutes)"
  ( cd exemplar-001 && ashlar \
      raw/exemplar-001-cycle-06.ome.tiff \
      raw/exemplar-001-cycle-07.ome.tiff \
      raw/exemplar-001-cycle-08.ome.tiff \
      --ffp illumination/exemplar-001-cycle-06-ffp.tif \
            illumination/exemplar-001-cycle-07-ffp.tif \
            illumination/exemplar-001-cycle-08-ffp.tif \
      --dfp illumination/exemplar-001-cycle-06-dfp.tif \
            illumination/exemplar-001-cycle-07-dfp.tif \
            illumination/exemplar-001-cycle-08-dfp.tif \
      --align-channel 0 \
      -o ../exemplar-001-registered.ome.tif )
else
  echo ">> exemplar-001-registered.ome.tif already present, skipping registration"
fi

if [ -d exemplar-001.zarr ]; then
  echo ">> exemplar-001.zarr already exists; remove it to rebuild"
else
  echo ">> Converting to SpatialData"
  python "$SCRIPT_DIR/ometiff_to_spatialdata.py" \
    --image exemplar-001-registered.ome.tif \
    --markers exemplar-001/markers.csv \
    --output exemplar-001.zarr \
    --image-key exemplar-001_image
fi

echo
echo ">> Done. Run the benchmark with:"
echo "   ZARR_PATH=$(pwd)/exemplar-001.zarr bash src/workflows/run_benchmark/test.sh"
