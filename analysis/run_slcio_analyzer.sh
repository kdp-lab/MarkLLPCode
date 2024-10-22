#!/bin/bash

# Ensure the script is called with at least three arguments.
if [ $# -lt 3 ]; then
  echo "Usage: $0 <number_of_chunks> <events_per_chunk> <sample_name1> [<sample_name2> ... <sample_nameN>]"
  exit 1
fi

# Get the first two arguments.
NUM_CHUNKS=$1
EVENTS_PER_CHUNK=$2

# Shift arguments so remaining ones are the sample names.
shift 2
SAMPLE_NAMES=("$@")  # Array of sample names.

# Loop over each sample name.
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
  echo "Processing sample: $SAMPLE_NAME"

  # Loop through each chunk for the current sample.
  for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
    echo "Running osg_slcio_analyzer.py with sample $SAMPLE_NAME, chunk $chunk"

    # Run the Python script with the current sample and chunk in the background.
    nohup python osg_slcio_analyzer.py --chunk $chunk \
      --eventsPerChunk $EVENTS_PER_CHUNK --sampleName $SAMPLE_NAME &

    # Check if the Python script started successfully.
    if [ $? -ne 0 ]; then
      echo "Error: Script failed to start for sample $SAMPLE_NAME, chunk $chunk"
      exit 1
    fi
  done

  # Wait for all background processes for this sample to finish.
  wait
  echo "All chunks processed for sample: $SAMPLE_NAME"

  # Create a list of all ROOT files for the current sample.
  ROOT_FILES=""
  for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
    ROOT_FILES+="${SAMPLE_NAME}_chunk${chunk}.root "
  done

  # Use hadd to merge all the ROOT files for the current sample.
  echo "Merging ROOT files into ${SAMPLE_NAME}.root"
  hadd "${SAMPLE_NAME}.root" $ROOT_FILES

  # Check if hadd ran successfully.
  if [ $? -ne 0 ]; then
    echo "Error: hadd failed to merge ROOT files for sample $SAMPLE_NAME"
    exit 1
  fi

  echo "Merged ROOT file created: ${SAMPLE_NAME}.root"
done

echo "All samples processed successfully!"
