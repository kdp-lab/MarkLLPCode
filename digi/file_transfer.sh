#!/bin/bash

#usage: nohup ./file_transfer.sh > transfer.log 2>&1 & 

# Set the source directory (current directory)
source_dir=$(pwd)

# Set the target directory (replace this with your desired target directory)
target_dir="/ospool/uc-shared/project/futurecolliders/larsonma/DigiMediumTiming10pBIB/"

# Run the loop indefinitely
while true; do
    # Find all files matching the pattern 'filename*' in the source directory
    for file in "$source_dir"/1000_1_digi*; do
        # Check if the file exists
        if [ -e "$file" ]; then
            echo "Transferring $file to $target_dir"
            mv "$file" "$target_dir"
        fi
    done

    # Sleep for a few seconds before checking again (adjust as needed)
    sleep 5
done
