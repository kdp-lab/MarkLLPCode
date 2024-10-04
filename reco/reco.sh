#!/bin/bash
echo $HOSTNAME
echo "<<<Singularity ENVIRONMENT:" $SINGULARITY_NAME
echo "<<<Setup some environment"
export MARLIN_DLL=$(realpath libMyBIBUtils.so):${MARLIN_DLL}
echo "source /opt/setup_mucoll.sh --> "
source /opt/setup_mucoll.sh
echo ">>>completed"
echo "<<<Check if we can find executables"
which ddsim
which Marlin
echo "<<<Check if input files were copied from the origin"
ls -lta 

# Initialize variables
input_file=""
output_directory=""
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

# Construct the nohup command
command="k4run reco_steer.py --LcioEvent.Files ${input_file}_timext_digi${proc_id}.slcio --outputFile ${input_file}_timext_reco${proc_id}.slcio --MatFile ${ACTS_MatFile} --TGeoFile ${ACTS_TGeoFile}"

# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command


echo "<<<copy that local file back to the origin"
echo "set stashcp client for non-OSG images"
# Copy finished sim file
cat  /etc/*-release  | grep VERSION_ID
export STASHCP=/cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/23/current/el8-x86_64/usr/bin/stashcp
$STASHCP -d ${input_file}_timext_reco${proc_id}.slcio osdf:///ospool/uc-shared/project/futurecolliders/larsonma/reco_osg_condor2/${input_file}_timext_reco${proc_id}.slcio
echo ">>> transfer completed"



echo "<<<Delete input files so they don't get transfered twice on exit"
rm -rf reco_steer.py
rm -rf ${input_file}_timext_digi${proc_id}.slcio
rm -rf ${input_file}_timext_reco${proc_id}.slcio
rm -rf libMyBIBUtils.so
echo ">>> Deletions complete. Test job complete"
