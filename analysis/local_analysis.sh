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

  # Determine the number of chunks per loop
  CHUNKS_PER_LOOP=20
  TOTAL_LOOPS=$((NUM_CHUNKS / CHUNKS_PER_LOOP))
  NUM_CORES=$(nproc)  # Total available CPU cores

  # Outer loop for each group of 20 chunks
  for (( loop=0; loop<TOTAL_LOOPS; loop++ )); do
    echo "Processing loop $((loop + 1)) of $TOTAL_LOOPS"

    # Inner loop for the chunks in this group
    for (( chunk=loop*CHUNKS_PER_LOOP; chunk<(loop+1)*CHUNKS_PER_LOOP; chunk++ )); do
      echo "Running osg_slcio_analyzer.py with sample $SAMPLE_NAME, chunk $chunk"

      # Assign the command to a specific core using `taskset`
      core=$((chunk % NUM_CORES))  # Rotate across available CPU cores
      nohup taskset -c $core python osg_slcio_analyzer.py --chunk $chunk \
        --eventsPerChunk $EVENTS_PER_CHUNK --sampleName $SAMPLE_NAME > ${SAMPLE_NAME}_${chunk}.txt &

      # Check if the Python script started successfully
      if [ $? -ne 0 ]; then
        echo "Error: Script failed to start for sample $SAMPLE_NAME, chunk $chunk"
        exit 1
      fi
    done

    # Optionally wait for processes in the current loop to complete
    wait
  done

echo "All loops completed successfully."

  # Wait for all background processes for this sample to finish.
  wait
  echo "All chunks processed for sample: $SAMPLE_NAME"

  # Create a list of all valid ROOT files for the current sample.
  ROOT_FILES=""
  for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
    ROOT_FILE="${SAMPLE_NAME}_chunk${chunk}.root"

    # Check if the ROOT file exists and is larger than 1 MB (1 MB = 1048576 bytes).
    if [ -f "$ROOT_FILE" ]; then
      FILE_SIZE=$(stat --format=%s "$ROOT_FILE")
      if [ "$FILE_SIZE" -ge 1048576 ]; then
        # Add the file to the list if it is valid.
        ROOT_FILES+="$ROOT_FILE "
      else
        echo "Warning: Skipping $ROOT_FILE (file is smaller than 1 MB)"
      fi
    else
      echo "Warning: Skipping $ROOT_FILE (file does not exist)"
    fi
  done

echo "All samples processed successfully!"
