#!/bin/bash
echo $HOSTNAME
echo "<<<Singularity ENVIRONMENT:" $SINGULARITY_NAME
echo "<<<Setup some environment"
echo "source source /opt/setup_mucoll.sh --> "
source /opt/setup_mucoll.sh
echo ">>>completed"
echo "<<<Check if we can find executables"
which ddsim
which Marlin
echo "<<<Check if input files were copied from the origin"
ls -lta 

# Initialize variables
input_file=""
chunks=""
n_events=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --inputFile)
            input_file="$2"
            shift 2
            ;;
        --chunks)
            chunks="$2"
            shift 2
            ;;
        --nEvents)
            n_events="$2"
            shift 2
            ;;
        --pId)
            proc_id="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Correct way to echo the proc_id
echo "Proc id: ${proc_id}"

# Correct variable access and calculation
events_per_chunk=$((n_events / chunks))
events_to_skip=$((proc_id * (n_events / chunks)))

# Corrected command assignment
command="ddsim --steeringFile steer_baseline.py --inputFile ${input_file}.hepmc --physics.pdgfile ${input_file}.tbl --outputFile ${input_file}_sim${proc_id}.slcio --numberOfEvents ${events_per_chunk} --skipNEvents ${events_to_skip}"

# Print or execute the constructed command
echo "Command: $command"
# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command


echo "<<<copy that local file back to the origin"
echo "set stashcp client for non-OSG images"
# Copy finished sim file
cat  /etc/*-release  | grep VERSION_ID
#export STASHCP=/cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/23/current/el8-x86_64/usr/bin/stashcp
#$STASHCP -d ${input_file}_sim${proc_id}.slcio osdf:///ospool/uc-shared/project/futurecolliders/larsonma/2500evSim/${input_file}_sim${proc_id}.slcio
echo ">>> transfer completed"

ls -lhtr

echo "<<<Delete input files so they don't get transferred twice on exit"
rm -rf steer_baseline.py
rm -rf $input_file.hepmc
rm -rf $input_file.tbl
#rm -rf ${input_file}_sim${proc_id}.slcio
echo ">>> Deletions complete. Test job complete"
