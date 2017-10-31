#!/usr/bin/env bash
# This script will run the steps for probtrackx in FSL's Diffusion Toolbox
# WITH PARALLEL JOBS
#
# ------------------------------------------
# ------------------STEPS-------------------
# ------------------------------------------
# 0) Create a subID.txt file in $script_dir with all subject IDs
# 1) Perform Eddy Correction and calculate tensors 
# 2) Create file structure for bedpostX
# 3) Run bedpostX
# 4) Run FLIRT to generate tranformation matrices (xfms)
#	a) 2mm (Default)
#	b) 1mm
# 5) Transform atlas and/or seed/target ROIs to use in probtrackX
#	a) FIRST segmentation
#	b) AAL Atlas
# 6) Run probtrackX:
# 	a) To generate a NxN connectivity matrix
#	b) To run a seed2target analysis
# 7) Normalize, threshold, binarize fdt_paths and extract mean FA from each tract
# ------------------------------------------
# ------------------------------------------
# ------------------------------------------


########run preprocessing in parallel########
FSLPARALLEL=1;  export FSLPARALLEL
MAXJOBS=12
sleeptime=0
function & waitforjobs {
	while [ $(jobs -p | wc -l) -ge $MAXJOBS ]; do
		echo "@$MAXJOBS jobs, sleeping $sleeptime s"
		jobs | sed 's/^/\t/'
		sleep $sleeptime
	done
}
########run preprocessing in parallel########


#=====================0=====================#
#===========CHANGE THESE OPTIONS============#
analysis_dir=/home/local/AD3/srhoads/w8_DTI_PGS/
script_dir="$analysis_dir"/scripts/
#=====================0=====================#


#=====================1=====================#
#========Convert raw dicoms to nii==========#
#======Also generates bval and bvecs========#
#=========Performs Eddy Correction==========#
#========BET on Eddy corrected image========#
#=============Calculate Tensors=============#
for i in $(cat "$script_dir"/all_DTIsubs.txt);
    do
	"$scriptdir"/preprocDTI1.bash "$i" &
	& waitforjobs
done
wait
#========Final files should include:========#
#	dti_SO.nii.gz							#
#	dti_FA.nii.gz							#
#	dti_MD.nii.gz							#
#	dti_L1-L3.nii.gz						#
#	dti_V1-V3.nii.gz						#
#=====================1=====================#
											

#=====================2=====================#
#											#
#		Will add soon!						#	
#	bval txt -> ~/analysis_dir/subID/bvals	#
#	bvec txt -> ~/analysis_dir/subID/bvecs	#
#	BETTED mprage -> ~/$analysis_dir/		#
#				../subID/BETmprage.nii.gz	#
#	eddy -> ~/$analysis_dir/data.nii.gz		#
#	bet_eddy_mask -> ~/$analysis_dir/		#
#				../nodif_brain_mask.nii.gz	#
#	bet_eddy.nii.gz -> ~/$analysis_dir/		#
#					../nodif.nii.gz			#
#											#
#	You can perform a check using 			#
#	'bedpostx_datacheck' to see if file 	#
#	structure is good for running bedpostX	#
#											#
#=====================2=====================#


#=====================3=====================#
#===============Run bedpostX================#
for j in $(cat "$script_dir"/all_DTIsubs.txt);
    do
    echo "$j"
	bedpostx $analysis_dir/"$j"/ -n 2 -model 1
	& waitforjobs
done
#=====================3=====================#

#=====================4=====================#
#=================Run FLIRT=================#
#=============Generates xmfs in=============#
#======$analysis_dir/*.bedpostX/xfms========#

# A) DEFAULT: generate xfms using 1mm standard space
echo "xfms using MNI152_T1_1mm_brain"
for k in $(cat "$script_dir"/all_DTIsubs.txt);
    do
    echo "$k"
    "$script_dir"/runSeed2Target_probtrackx5000/runFLIRT.sh "$k" & 
    waitforjobs
done

# B) Can also change to "MNI152_T1_2mm_brain if you want to generate xfms for 2mm standard space
#			for use with FIRST subcortical segmentation
#=====================4=====================#


#=====================5a====================#
#============FIRST segmentation=============#
#===puts subcortical ROIs into diff space===#
#============uses 1mm MNI space=============#
#============Binarize ROI masks=============#
for l in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$k"
    "$script_dir"/runSeed2Target_probtrackx5000/runFLIRT.sh "$k" &
	waitforjobs
done
#=====================5a====================#

#=====================5b====================#
#==========AAL atlas Parcellation===========#
#											#
#	Adding later							#
#											#
#=====================5b====================#

#=====================6a====================#
#======run probtrackX: NxN Connectivity=====#
#											#
#	Adding later							#
#											#
#=====================6a====================#

#=====================6b====================#
#=======run probtrackX: seed2target(s)======#
# Make target mask list in MATLAB:
matlab makeTargetList.m

#=====================6b====================#

#=====================6c====================#
#========run probtackX: N Seed Masks========#
# Make seed mask list in MATLAB:
matlab makeSeekMask.m

for nn in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	"$script_dir"/runSeed2Target_probtrackx5000/runProbtrackX_2MaskSeeding.sh "$nn" &
	waitforjobs
done
#=====================6c====================#

#=====================7=====================#
#=========Normalize fdt_paths.nii.gz========#
#=========Threshold fdt_paths.nii.gz========#
#=========Binarize fdt_paths.nii.gz=========#
#====Extract mean FA and mean Intensity=====#
for subID in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	cd "$analysis_dir"/"$subID".probtrackx5000/MPFC2Striatum/

	# calculate max intensity to create percentage map per subject
	fslstats fdt_paths.nii.gz -R >> intensityRange.txt
	maxIntensity=$(awk '{print $2}' intensityRange.txt)

	# generate normalized fdt_paths 
	fslmaths fdt_paths.nii.gz -div "$maxIntensity" norm_fdt_paths.nii.gz

	# generate thresholded mask at 50%
	fslmaths norm_fdt_paths.nii.gz -thr 0.5 thr50_norm_fdt_paths.nii.gz

	# binarize into mask
	fslmaths thr50_norm_fdt_paths.nii.gz -bin bin_thr50_norm_fdt_paths.nii.gz

	#calculate mean intensity in normalized tract
	fslstats thr50_norm_fdt_paths.nii.gz -M >> meanIntensity.txt

	# calculate mean FA within binarized mask
	fslstats /mnt/datashare/Data/Raj/dti_8/"$subID"/"$subID"_FA.nii.gz -k bin_thr50_norm_fdt_paths.nii.gz -M >> meanFA.txt

	echo "Thresholded at 50%"
	echo "$subID",$(cat "meanFA.txt"),$(cat "meanIntensity.txt")

	# generate thresholded mask at 1%
	fslmaths norm_fdt_paths.nii.gz -thr 0.01 thr01_norm_fdt_paths.nii.gz

	# binarize into mask
	fslmaths thr01_norm_fdt_paths.nii.gz -bin bin_thr01_norm_fdt_paths.nii.gz

	#calculate mean intensity in normalized tract
	fslstats thr01_norm_fdt_paths.nii.gz -M >> meanIntensity_thr01.txt

	# calculate mean FA within binarized mask
	fslstats /mnt/datashare/Data/Raj/dti_8/"$subID"/"$subID"_FA.nii.gz -k bin_thr01_norm_fdt_paths.nii.gz -M >> meanFA_thr01.txt

	echo "Thresholded at 1%"
	echo "$subID",$(cat "meanFA_thr01.txt"),$(cat "meanIntensity_thr01.txt")
done
#=====================7=====================#
