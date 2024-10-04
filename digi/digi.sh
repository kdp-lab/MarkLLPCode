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

#mkdir sim_mp_pruned
#mkdir sim_mm_pruned

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

export STASHCP=/cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/23/current/el8-x86_64/usr/bin/stashcp

# Construct the path to the file you want to copy based on random_number
#file_to_copy="osdf:///ospool/uc-shared/project/muoncollider/tutorial2024/MuColl_v1/BIB10TeV/sim_mp_pruned/"

# Use stashcp to copy the file to the worker's local directory
#$STASHCP -d -r $file_to_copy /sim_mp_pruned

# Construct the path to the file you want to copy based on random_number
#file_to_copy="osdf:///ospool/uc-shared/project/muoncollider/tutorial2024/MuColl_v1/BIB10TeV/sim_mm_pruned/"

# Use stashcp to copy the file to the worker's local directory
#$STASHCP -d -r $file_to_copy /sim_mm_pruned

ls sim_mm_pruned

ls sim_mp_pruned

: << 'COMMENT'
### TRANSFER NECESSARY NUMBER OF BIB FILES
# Set N (number of random numbers) and n (upper limit)
N=${bib}
n=768



# Initialize an array to store unique random numbers
unique_numbers=()

# Loop through N random numbers
for ((i=0; i<N; i++))
do
    while true
    do
        # Generate a random number between 0 and n
        random_number=$((RANDOM % n))

        # Check if the number is already in the array
        if [[ ! " ${unique_numbers[@]} " =~ " ${random_number} " ]]; then
            # If not, store it in the array and break out of the while loop
            unique_numbers+=("$random_number")
            echo "Random number $i: $random_number"

            # Construct the path to the file you want to copy based on random_number
            file_to_copy="osdf:///ospool/uc-shared/project/muoncollider/tutorial2024/MuColl_v1/BIB10TeV/sim_mp_pruned/BIB_sim_${random_number}.slcio"

            # Use stashcp to copy the file to the worker's local directory
            $STASHCP -d $file_to_copy /sim_mp_pruned/file_${random_number}.ext

            # Construct the path to the file you want to copy based on random_number
            file_to_copy="osdf:////ospool/uc-shared/project/muoncollider/tutorial2024/MuColl_v1/BIB10TeV/sim_mm_pruned/BIB_sim_${random_number}.slcio"

            # Use stashcp to copy the file to the worker's local directory
            $STASHCP -d $file_to_copy /sim_mm_pruned/BIB_sim_${random_number}.slcio

            echo "File for random number $random_number copied."

            break
        fi
    done
done
COMMENT
# Construct the nohup command

MUPLUS="sim_mp_pruned/"
MUMINUS="sim_mm_pruned/"

command="k4run digi_steer.py --LcioEvent.Files ${input_file}_sim${proc_id}.slcio --outputFile ${input_file}_digi${proc_id}.slcio --doOverlayFull --OverlayFullPathToMuPlus ${MUPLUS} --OverlayFullPathToMuMinus ${MUMINUS} --OverlayFullNumberBackground ${bib}"
# To not run with BIB set bib = 0
# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command

echo "<<<copy that local file back to the origin"
echo "set stashcp client for non-OSG images"
# Copy finished sim file
cat  /etc/*-release  | grep VERSION_ID

$STASHCP -d ${input_file}_digi${proc_id}.slcio osdf:///ospool/uc-shared/project/futurecolliders/larsonma/digi_osg_condor/${input_file}_digi${proc_id}.slcio
echo ">>> transfer completed"



echo "<<<Delete input files so they don't get transfered twice on exit"
rm -rf digi_steer.py
rm -rf ${input_file}_sim${proc_id}.slcio
rm -rf ${input_file}_digi${proc_id}.slcio
rm -rf sim_mp_pruned
rm -rf sim_mm_pruned
echo ">>> Deletions complete. Test job complete"
