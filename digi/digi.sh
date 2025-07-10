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

wget https://dl.pelicanplatform.org/7.15.3/pelican_Linux_x86_64.tar.gz
gunzip pelican_Linux_x86_64.tar.gz 
tar -xvmf pelican_Linux_x86_64.tar 
export PATH=${PWD}/pelican-7.15.3:$PATH

target_federation=osg-htc.org
origin_directory=/ospool/uc-shared/project/futurecolliders


# Initialize variables
input_file=""
output_directory=""
chunks=""
n_events=""
bib=""

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
        --BIB)
            bib="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

echo "Proc id: ${proc_id}"

# Construct the nohup command

MUPLUS="/cvmfs/public-uc.osgstorage.org/ospool/uc-shared/public/futurecolliders/MuColl_v1/BIB10TeV/sim_mp_pruned/"
MUMINUS="/cvmfs/public-uc.osgstorage.org/ospool/uc-shared/public/futurecolliders/MuColl_v1/BIB10TeV/sim_mm_pruned/"

command="k4run digi_steer.py --LcioEvent.Files ${input_file}_sim${proc_id}.slcio --outputFile ${input_file}_digi${proc_id}.slcio --doOverlayFull --OverlayFullPathToMuPlus ${MUPLUS} --OverlayFullPathToMuMinus ${MUMINUS} --OverlayFullNumberBackground ${bib}"
# To not run with BIB set bib = 0
# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command

export STASHCP=/cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/23/current/el8-x86_64/usr/bin/stashcp

echo "<<<copy that local file back to the origin"
echo "set stashcp client for non-OSG images"
# Copy finished sim file
cat  /etc/*-release  | grep VERSION_ID

pelican object put -d ${input_file}_digi${proc_id}.slcio pelican://${target_federation}${origin_directory}/larsonma/DigiTestJobs/${input_file}_digi${proc_id}.slcio
#$STASHCP -d ${input_file}_digi${proc_id}.slcio osdf:///ospool/uc-shared/project/futurecolliders/larsonma/DigiMediumTiming10pBIB/${input_file}_digi${proc_id}.slcio
#transfer_output_remaps = "${input_file}_digi${proc_id}.slcio=osdf:///ospool/uc-shared/project/futurecolliders/larsonma/DigiMediumTiming10pBIB/${input_file}_digi${proc_id}.slcio"
echo ">>> transfer completed"

echo ""
ls -l

echo "<<<Delete input files so they don't get transfered twice on exit"
rm -rf digi_steer.py
rm -rf ${input_file}_sim${proc_id}.slcio
rm -rf ${input_file}_digi${proc_id}.slcio
echo ">>> Deletions complete. Test job complete"
