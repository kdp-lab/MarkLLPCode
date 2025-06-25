#!/bin/bash
echo $HOSTNAME
echo "<<<Singularity ENVIRONMENT:" $SINGULARITY_NAME
echo "<<<Setup some environment"
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
command="python analysis.py --chunk ${proc_id} --sampleName ${input_file}"

# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command


echo "<<<copy that local file back to the origin"
echo "set stashcp client for non-OSG images"
# Copy finished sim file
cat  /etc/*-release  | grep VERSION_ID
export STASHCP=/cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/23/current/el8-x86_64/usr/bin/stashcp
$STASHCP -d ${input_file}_analysis${proc_id}.root /scratch/larsonma/tutorial2024/LLPStudies/MarkLLPCode/analysis/AnalysisJobsOutputNoBib/${input_file}_analysis${proc_id}.root
echo ">>> transfer completed"



echo "<<<Delete input files so they don't get transfered twice on exit"
rm -rf analysis.py
rm -rf ${input_file}_reco${proc_id}.slcio
echo ">>> Deletions complete. Test job complete"
