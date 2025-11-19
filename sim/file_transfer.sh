#!/bin/bash

#usage: nohup ./file_transfer.sh > transfer.log 2>&1 &
#get the source directory (current directory)
source_dir=$(pwd)
# Set the target directory (replace this with your desired target directory)
target_dir="/ospool/uc-shared/project/futurecolliders/larsonma/SimTest"

#run loop indefinitely
while true; do
	for file in "$source_dir"/4000_10_sim*; do
		[ -e "$file" ] || continue
		echo "Moving $file -> $target_dir"
		mv -- "$file" "$target_dir"
	done	

	sleep 5

done




