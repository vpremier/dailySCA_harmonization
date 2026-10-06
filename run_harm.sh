#!/bin/bash

# Exit if any command fails
set -e

# Initialize micromamba shell (only needed once per session/script)
eval "$(micromamba shell hook --shell bash)"

# Activate environment
echo "Activating micromamba environment 'snowmap'..."
micromamba activate snowmap

# Path to your Python script and config
SCRIPT_PATH="./main_correction.py"
CONFIG_PATH="./config.json"

for year in 2021 1314 1415 1516 1617 1718 1819 1920 2122 2223; do
  new_val="hy${year}"
  for catch in Area06; do
  # for catch in Area01 Area02 Area03 Area04 Area05 Area06 Area07 Area08 Area09 Area10; do
  
    echo "=============================================="
    echo "Processing $new_val - $catch"
    echo "=============================================="

    # Update JSON fields dynamically
    jq --arg v "$new_val" --arg c "$catch" '
      .hy_xxxx = $v
      | .catchment = $c
      | .dirname   |= sub("Area[0-9]+"; $c)
      | .outdir    |= sub("Area[0-9]+"; $c)
      | .DEM_path  |= sub("Area[0-9]+"; $c)
      | .temp_dir  |= sub("Area[0-9]+"; $c)
      | .era5_dir  |= sub("Area[0-9]+"; $c)
    ' "$CONFIG_PATH" > tmp.json && mv tmp.json "$CONFIG_PATH"

    echo "Updated config.json:"
    echo "  hy_xxxx = $new_val"
    echo "  catchment = $catch"

    # Run the Python script
    echo "Running: python $SCRIPT_PATH $CONFIG_PATH"
    python "$SCRIPT_PATH" "$CONFIG_PATH"

  done
done

echo "✅ All processing completed successfully."

