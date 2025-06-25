
# Get the first two arguments.
NUM_CHUNKS=500

# Shift arguments so remaining ones are the sample names.
SAMPLE_NAME=4000_10  # Array of sample names.


# Create a list of all valid ROOT files for the current sample.
ROOT_FILES=""
for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
  ROOT_FILE="${SAMPLE_NAME}_analysis${chunk}.root"

  # Check if the ROOT file exists and is larger than 1 MB (1 MB = 1048576 bytes).
  if [ -f "$ROOT_FILE" ]; then
    FILE_SIZE=$(stat --format=%s "$ROOT_FILE")
    if [ "$FILE_SIZE" -lt 1024 ]; then
      echo "Warning: Skipping $ROOT_FILE (file is smaller than 1 KB)"
      else
        # Add the file to the list if it is valid.
        ROOT_FILES+="$ROOT_FILE "
    fi
  else
    echo "Warning: Skipping $ROOT_FILE (file does not exist)"
  fi
done

# Check if there are any valid ROOT files to merge.
if [ -z "$ROOT_FILES" ]; then
  echo "Error: No valid ROOT files found for sample $SAMPLE_NAME"
  exit 1
fi

# Use hadd to merge all the valid ROOT files for the current sample.
echo "Merging ROOT files into ${SAMPLE_NAME}.root"
hadd "${SAMPLE_NAME}_reco.root" $ROOT_FILES

mv ${SAMPLE_NAME}_reco.root openHouseRootFiles/${SAMPLE_NAME}_medium_velores_reco.root

#rm -rf $ROOT_FILES

echo "All samples processed successfully!"
