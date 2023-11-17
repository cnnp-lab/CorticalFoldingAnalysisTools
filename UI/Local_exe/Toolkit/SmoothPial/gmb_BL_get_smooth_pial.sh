#!/bin/bash

# This script gets the smooth pial files (runs the first four steps from the localGI command in Freesurfer). 
# Requires a tab delimited txt file that lists the subject IDs in a row (subjects.txt). 
# Require FreeSurfer and Matlab to be installed. 

# manually load and source FreeSurfer 
export FREESURFER_HOME="$(cat fs_path.txt)"
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Define the Matlab directory
export matlab="$(cat mtl_path.txt)"

# Get the list of subjects to process (listed in txt file).
subs="$(cat subjects.txt)"

# for each subject...
#for subj_id in $subs; do 

# ...for each hemisphere...
	for hemisphere in lh rh; do

		#...run first four steps from lGI process to get pial-outer-smoothed surface.
		mris_fill -c -r 1 ${subs}/surf/${hemisphere}.pial ${subs}/surf/${hemisphere}.pial.filled.mgz
		
		${matlab} -nodisplay -nosplash -nodesktop -r "make_outer_surface('${subs}/surf/${hemisphere}.pial.filled.mgz', 15, '${subs}/surf/${hemisphere}.pial-outer'); exit;"
		
		mris_extract_main_component ${subs}/surf/${hemisphere}.pial-outer ${subs}/surf/${hemisphere}.pial-outer-main
		
		mris_smooth -nw -n 30 ${subs}/surf/${hemisphere}.pial-outer-main ${subs}/surf/${hemisphere}.pial-outer-smoothed

	done
done


