#!/bin/bash

# This script gets the smooth pial files (runs the first four steps from the localGI command in Freesurfer). 
# Requires a tab delimited txt file that lists the subject IDs in a row (subjects.txt). 
# Require FreeSurfer and Matlab to be installed. 

# manually load and source FreeSurfer 
export FREESURFER_HOME=/Applications/freesurfer/7.4.1
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# set subject directory to the the folder where recon-all output is stored
export SUBJECTS_DIR=/Users/guillermo/Sync/Data/04_CorticalFolding_Interface/gmb_Trial/FreeSurfer
#echo $SUBJECTS_DIR 

# Get the list of subjects to process (listed in txt file).
export SUBJ_ID=/Users/guillermo/Sync/Projects/04_CorticalFolding_Interface/gmb_CotFold_Local_V2/Local_exe/Toolkit/SmoothPial
subs="$(cat ${SUBJ_ID}/subjects.txt)"
#echo $subs

# for each subject...
for subj_id in $subs; do 

	# ...for each hemisphere...
	for hemisphere in lh rh; do

		#...run first four steps from lGI process to get pial-outer-smoothed surface.
		mris_fill -c -r 1 ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial.filled.mgz
		
		matlab -nodisplay -nosplash -nodesktop -r "make_outer_surface('${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial.filled.mgz', 15, '${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial-outer'); exit;"
		
		#mris_extract_main_component ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial-outer ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial-outer-main
		
		#mris_smooth -nw -n 30 ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial-outer-main ${SUBJECTS_DIR}/${subj_id}/surf/${hemisphere}.pial-outer-smoothed

	done
done
