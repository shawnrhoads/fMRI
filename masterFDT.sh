#!/usr/bin/env bash
# This script will run the steps for probtrackx in FSL's Diffusion Toolbox
# WITHOUT PARALLEL JOBS --- to run in parallel, refer to masterFDT_parallel.bash
#
# ------------------------------------------
# ------------------STEPS-------------------
# ------------------------------------------
# 0) Create a subID.txt file in $script_dir with all subject IDs
	# a) "$analysis_dir"/ should contain directories with subject IDs and raw dicom/nifti files from DTI
# 1) Perform Eddy Correction and calculate tensors 
# 2) Create file structure for bedpostX
# 3) Run bedpostX
# 4) Run FLIRT to generate tranformation matrices (xfms)
#	a) MNI152_T1_1mm (Default) but can also use MNI152_T1_2mm
# 5) Transform atlas and/or seed/target ROIs to use in probtrackX
#	a) AAL Atlas
#	b) FIRST segmentation
#	c) ROI sphere seed to diffusion space
# 6) Run probtrackX:
# 	a) To generate a NxN connectivity matrix
#	b) To run a seed2target analysis
# 7) Normalize, threshold, binarize fdt_paths and extract mean FA from each tract
# ------------------------------------------
# ------------------------------------------
# ------------------------------------------


#=====================0=====================#
#===========CHANGE THESE OPTIONS============#
analysis_dir=/home/local/AD3/srhoads/w8_DTI_PGS
script_dir="$analysis_dir"/scripts
#=====================0=====================#


#=====================1=====================#
#========Convert raw dicoms to nii==========#
#======Also generates bval and bvecs========#
#=========Performs Eddy Correction==========#
#========BET on Eddy corrected image========#
#=============Calculate Tensors=============#
for i in $(cat "$script_dir"/all_DTIsubs.txt);
    do
    # Convert DTI dicoms to nii (if relevant)
	dcm2nii "$analysis_dir"/"$i"/raw/

	# I also wrote a script in MATLAB using SPM to move and convert dicoms to nii
	## Move and convert mprage.img files to mprage.nii
	## wavePath = '/mnt/datashare/Data/UPitt/wave8/';
	## subID=dlmread('/home/local/AD3/srhoads/w8_DTI_PGS/scripts/all_DTIsubs.txt');
	## anatPath= '/Anat/raw/';
	## dtiPath='/mnt/datashare/Data/Shawn/w8_PGSE_DTI/';
	## for i = 1:length(subID)
	##     try 
	##         %---convert .img to .nii using SPM---%
	##         f=strcat(wavePath,subID{i},anatPath,'s',subID{i},'*.img');
	##         file=dir(f);
	##         filepath=strcat(wavePath,subID{i},anatPath,file.name);
	##         for j=1:size(filepath,1)
	##             input = deblank(filepath(j,:));
	##             [pathstr,fname,ext] = fileparts(input);
	##             output = strcat(dtiPath,subID{i},'/mprage.nii'); %put in DTI directory
	##             V=spm_vol(input);
	##             ima=spm_read_vols(V);
	##             V.fname=output;
	##             spm_write_vol(V,ima);
	##         end
	##         %---convert .img to .nii using SPM---%
	##     catch err
	##         print(err)
	##     end        
	## end


	#Archive dicoms
	nifti="$analysis_dir"/"$i"/raw/*.nii.gz
	chmod -R 777 "$analysis_dir"/"$i"/raw/*

	#Eddy Correction
	eddy_correct "$nifti" "$analysis_dir"/"$i"/data --0

	#BET on Eddy Corrected Image (nodif.nii.gz) + Mask (nodif_mask.nii.gz)
	bet "$analysis_dir"/"$i"/data.nii.gz "$analysis_dir"/"$i"/nodif -m -f 0.2

	#Calculate tensors
	dtifit -k "$analysis_dir"/"$i"/data.nii.gz -o "$analysis_dir"/"$i"/dti -m "$analysis_dir"/"$i"/nodif_mask.nii.gz -r "$analysis_dir"/"$i"/raw/*.bvec -b "$analysis_dir"/"$i"/raw/*.bval
done
#===========================================#
# Final files should include:
	# <basename>_V1 - 1st eigenvector
	# <basename>_V2 - 2nd eigenvector
	# <basename>_V3 - 3rd eigenvector
	# <basename>_L1 - 1st eigenvalue
	# <basename>_L2 - 2nd eigenvalue
	# <basename>_L3 - 3rd eigenvalue
	# <basename>_MD - mean diffusivity
	# <basename>_FA - fractional anisotropy
	# <basename>_MO - mode of the anisotropy (oblate ~ -1; isotropic ~ 0; prolate ~ 1)
	# <basename>_S0 - raw T2 signal with no diffusion weightingdti_SO.nii.gz
#=====================1=====================#
											

#=====================2=====================#
#======Make file structure for bedpostX=====#
#	bval txt -> ~/analysis_dir/subID/bvals	#
#	bvec txt -> ~/analysis_dir/subID/bvecs	#
#	BETTED mprage -> ~/$analysis_dir/		#
#				../subID/BETmprage.nii.gz	#
#	eddy -> ~/$analysis_dir/data.nii.gz		#
#	bet_eddy_mask -> ~/$analysis_dir/		#
#				../nodif_brain_mask.nii.gz	#
#	bet_eddy.nii.gz -> ~/$analysis_dir/		#
#					../nodif.nii.gz			#
#===========================================#
# Move bvals and bvecs up one directory without extension
for ii in $(cat "$script_dir"/all_DTIsubs.txt);
    do
	cp "$analysis_dir"/"$ii"/raw/*.bval "$analysis_dir"/"$ii"/bvals;
	cp "$analysis_dir"/"$ii"/raw/*.bvec "$analysis_dir"/"$ii"/bvecs;
done
#===========================================#
#											#
#	You can perform a check using 			#
#	'bedpostx_datacheck' to see if file 	#
#	structure is good for running bedpostX 	#
#	bedpostX expects to find:				#
#		- bvals								#
#		- bvecs								#
#		- data 								#
#		- nodif_brain_mask					#
#			in subject directory			#
#				"$analysis_dir"/"$i"/		#
#											#
#=====================2=====================#


#=====================3=====================#
#===============Run bedpostX================#
# Consider running this step in PARALLEL, see masterFDT_parallel.bash
for j in $(cat "$script_dir"/all_DTIsubs.txt);
    do
    echo "$j"
	bedpostx $analysis_dir/"$j"/ -n 2 -model 1
done
#=====================3=====================#


#=====================4=====================#
#=================Run FLIRT=================#
#=============Generates xmfs in=============#
#======$analysis_dir/*.bedpostX/xfms========#

# A) DEFAULT: generate xfms using 1mm standard space (for use with FIRST subcortical segmentation)
echo "xfms using MNI152_T1_1mm_brain"
for k in $(cat "$script_dir"/all_DTIsubs.txt);
    do
    echo "$k"
	#IMPORTANT NOTE: need to use BETTED MPRAGE brain (MNI152_T1_1mm_brain.nii.gz)
	echo "betting"
	bet2 "$analysis_dir"/"$k"/mprage.nii "$analysis_dir"/"$k"/BETmprage.nii
	echo "diff2str"
	flirt -in "$analysis_dir"/"$k".bedpostX/nodif_brain.nii.gz -ref "$analysis_dir"/"$k"/BETmprage.nii.gz -omat "$analysis_dir"/"$k".bedpostX/xfms/diff2str.mat 
	echo "str2diff"
	convert_xfm -omat "$analysis_dir"/"$k".bedpostX/xfms/str2diff.mat -inverse "$analysis_dir"/"$k".bedpostX/xfms/diff2str.mat
	echo "str2standard"
	flirt -in "$analysis_dir"/"$k"/mprage.nii -ref /usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz -omat "$analysis_dir"/"$k".bedpostX/xfms/str2standard.mat
	echo "standard2str"
	convert_xfm -omat "$analysis_dir"/"$k".bedpostX/xfms/standard2str.mat -inverse "$analysis_dir"/"$k".bedpostX/xfms/str2standard.mat
	#IMPORTANT NOTE: need to use BETTED MNI brain (MNI152_T1_1mm_brain.nii.gz)
	echo "diff2standard"
	flirt -in "$analysis_dir"/"$k".bedpostX/nodif_brain.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -omat "$analysis_dir"/"$k".bedpostX/xfms/diff2standard.mat 
	echo "standard2diff"
	convert_xfm -omat "$analysis_dir"/"$k".bedpostX/xfms/standard2diff.mat -inverse "$analysis_dir"/"$k".bedpostX/xfms/diff2standard.mat
	echo "Done"
done

# B) Can also change to "MNI152_T1_2mm_brain if you want to generate xfms for 2mm standard space

#=====================4=====================#


#=====================5a====================#
#==========AAL atlas Parcellation===========#
# Works with:
#	- AAL atlas in "$analysis_dir"/AAL_masks/
#	- Oriented AAL as "$analysis_dir"/AAL_masks/aal2_orient.nii.gz
#	- AAL region list in "$script_dir"/aal_list.txt

for l in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$l"
	nodif="$analysis_dir"/"$l".bedpostX/nodif_brain.nii.gz
	aal2="$analysis_dir"/AAL_masks/aal2_orient.nii.gz
	mkdir "$analysis_dir"/"$l".bedpostX/ROIs/AAL_transform/
	transformdir="$analysis_dir"/"$l".bedpostX/ROIs/AAL_transform/

	# Generate aal2diff and diff2aal xfms
	flirt -in $nodif -ref $aal2 "$analysis_dir"/"$l".bedpostX/xfms/diff2aal.mat 
	convert_xfm -omat "$analysis_dir"/"$l".bedpostX/xfms/aal2diff.mat -inverse "$analysis_dir"/"$l".bedpostX/xfms/diff2aal.mat

	for region in $(cat "$script_dir"/aal_list.txt)
		do
		applywarp --ref=$nodif --in="$analysis_dir"/AAL_masks/"$region"_AALmask.nii.gz --postmat="$analysis_dir"/"$l".bedpostX/xfms/diff2aal.mat --out="$transformdir"/"$region"_AALmask-diff
		fslmaths "$transformdir"/"$region"_AALmask-diff.nii.gz -bin "$transformdir"/bin_"$region"_AALmask-diff.nii.gz
	done
done
#=====================5a====================#


#=====================5b====================#
#============FIRST segmentation=============#
#===puts subcortical ROIs into diff space===#
#============uses 1mm MNI space=============#
#============Binarize ROI masks=============#
for ll in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$ll"
	mkdir "$analysis_dir"/"$ll".bedpostX/ROIs/
	echo "running FIRST"
	# Segment Left Caudate
	run_first -i "$analysis_dir"/"$ll".bedpostX/nodif_brain.nii.gz -t "$analysis_dir"/"$ll".bedpostX/xfms/diff2standard.mat -n 30 -o "$analysis_dir"/"$ll".bedpostX/ROIs/L_caudate_n30 -m /usr/local/fsl/data/first/models_336_bin/L_Caud_bin.bmv
	# Segment Right Caudate
	run_first -i "$analysis_dir"/"$ll".bedpostX/nodif_brain.nii.gz -t "$analysis_dir"/"$ll".bedpostX/xfms/diff2standard.mat -n 30 -o "$analysis_dir"/"$ll".bedpostX/ROIs/R_caudate_n30 -m /usr/local/fsl/data/first/models_336_bin/R_Caud_bin.bmv
	# Segment Left Accumbens
	run_first -i "$analysis_dir"/"$ll".bedpostX/nodif_brain.nii.gz -t "$analysis_dir"/"$ll".bedpostX/xfms/diff2standard.mat -n 50 -o "$analysis_dir"/"$ll".bedpostX/ROIs/L_accumbens_n50 -m /usr/local/fsl/data/first/models_336_bin/L_Accu_bin.bmv
	# Segment Right Accumbens
	run_first -i "$analysis_dir"/"$ll".bedpostX/nodif_brain.nii.gz -t "$analysis_dir"/"$ll".bedpostX/xfms/diff2standard.mat -n 50 -o "$analysis_dir"/"$ll".bedpostX/ROIs/R_accumbens_n50 -m /usr/local/fsl/data/first/models_336_bin/R_Accu_bin.bmv

	echo "combining ROIs"
	# R_acc + R_cau = R_striatum
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/R_accumbens_n50.nii.gz -add "$analysis_dir"/"$ll".bedpostX/ROIs/R_caudate_n30.nii.gz "$analysis_dir"/"$ll".bedpostX/ROIs/R_striatum.nii.gz
	# L_acc + L_cau = L_striatum
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/L_accumbens_n50.nii.gz -add "$analysis_dir"/"$ll".bedpostX/ROIs/L_caudate_n30.nii.gz "$analysis_dir"/"$ll".bedpostX/ROIs/L_striatum.nii.gz
	# R_striatum + L_striatum = bilateral_striatum
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/R_striatum.nii.gz -add "$analysis_dir"/"$ll".bedpostX/ROIs/L_striatum.nii.gz "$analysis_dir"/"$ll".bedpostX/ROIs/bilateral_striatum.nii.gz

	# Binarize ROIs into masks
	echo "binarizing ROIs"
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/bilateral_striatum.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/bilateral_striatum_bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/L_striatum.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/L_striatum_bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/R_striatum.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/R_striatum_bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/L_accumbens_n50.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/L_accumbens_n50bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/R_accumbens_n50.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/R_accumbens_n50bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/R_caudate_n30.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/R_caudate_n30bin.nii.gz
	fslmaths "$analysis_dir"/"$ll".bedpostX/ROIs/L_caudate_n30.nii.gz -bin "$analysis_dir"/"$ll".bedpostX/ROIs/L_caudate_n30bin.nii.gz
done
#=====================5b====================#


#=====================5c====================#
#=====Put ROI sphere in diffusion space=====#
# From MNI_1mm standard to diff:
for lll in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$lll"
	flirt -in ~/MNI1mm_sphere_10mm-3_53_3.nii.gz -ref "$analysis_dir"/"$lll".bedpostX/nodif_brain.nii.gz -applyxfm -init "$analysis_dir"/"$lll".bedpostX/xfms/standard2diff.mat -out "$analysis_dir"/"$lll".bedpostX/ROIs/MPFCsphere-3_53_3-diff.nii

	flirt -in ~/MNI1mm_sphere_6mm-10-52-2.nii.gz -ref "$analysis_dir"/"$lll".bedpostX/nodif_brain.nii.gz -applyxfm -init "$analysis_dir"/"$lll".bedpostX/xfms/standard2diff.mat -out "$analysis_dir"/"$lll".bedpostX/ROIs/MPFCsphere-10_52_2-diff.nii
done
#=====================5c====================#


#=====================6a====================#
#======run probtrackX: NxN Connectivity=====#
# Consider running this step in PARALLEL, see masterFDT_parallel.bash

# Make mask_list in MATLAB
# matlab maskMaskList.m
# puts .txt file in *.bedpostX/ directory

probtrackx_type=matrix;

for m in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$m"
	probtrackx2 --samples="$analysis_dir"/"$m".bedpostX/merged --mask="$analysis_dir"/"$m".bedpostX/nodif_brain_mask --seed="$analysis_dir"/"$m".bedpostX/mask_list.txt --loopcheck --forcedir --network --omatrix1 -P 5000 -V 0 --dir="$analysis_dir"/"$m".probtrackX/"$probtrackx_type"
done
#=====================6a====================#


#=====================6b====================#
#=======run probtrackX: seed2target(s)======#
# Consider running this step in PARALLEL, see masterFDT_parallel.bash

# Make target mask list in MATLAB:
# matlab makeTargetList.m
# puts .txt file in *.bedpostX/ directory

probtrackx_type=seed2target;

for mm in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	echo "$mm"
	probtrackx2 -x "$analysis_dir"/"$mm".bedpostX/MPFCsphere_6mm-str.nii.gz  -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm="$analysis_dir"/"$mm".bedpostX/xfms/standard2diff.mat --forcedir --opd -s "$analysis_dir"/"$mm".bedpostX/merged -m "$analysis_dir"/"$mm".bedpostX/nodif_brain_mask  --dir="$analysis_dir"/"$mm".probtrackX/"$probtrackx_type"
done
#=====================6b====================#


#=====================6c====================#
#========run probtackX: N Seed Masks========#
# Consider running this step in PARALLEL, see masterFDT_parallel.bash

# Make seed mask list in MATLAB:
# matlab makeSeedMask.m
# puts .txt file in *.bedpostX/ directory

probtrackx_type=n_seed;

for nnn in $(cat "$script_dir"/all_DTIsubs.txt);
	do
	probtrackx2 --network -x "$analysis_dir"/"$nnn".bedpostX/MPFC2StriatumMaskList.txt -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s "$analysis_dir"/"$nnn".bedpostX/merged -m "$analysis_dir"/"$nnn".bedpostX/nodif_brain_mask --dir="$analysis_dir"/"$nnn".probtrackx5000/"$probtrackx_type"
done
#=====================6c====================#


#=====================7a====================#
#======FOR USE WITH PROBTRACKX: MATRIX======#

# Normalize all network matrices in MATLAB
	## for k = 1:length(subs)
	##     subID=subs{k};
	##     cd(strcat('/mnt/datashare/Data/Shawn/PGSE/wave8/',subID,'.probtrackx5000/matrix'))
	##     waytotal=load('waytotal');
	##     fdt_network_matrix=load('fdt_network_matrix');
	##     for i = 1:length(waytotal)
	##         for j = 1:length(fdt_network_matrix)
	##             norm_fdt_network_matrix(i,j) = fdt_network_matrix(i,j)/waytotal(i);
	##         end
	##     end
	##     save(strcat('norm_fdt_network_matrix.mat'), 'norm_fdt_network_matrix')
	##     clear waytotal fdt_network_matrix subID norm_fdt_network_matrix
	## end

# Put concat all normalized matrices in MATLAB
	## conn_mat = struct('',[]); %create struct to put everything in
	## subs={}; %list of subject IDs
	## for k = 1:length(subs)
	## 	cd(strcat('/mnt/datashare/Data/Shawn/PGSE/wave8/',subID,'.probtrackx5000/matrix'))
	##     load(strcat(subs{k},'_matrix_norm.mat'))
	##     fieldname=strcat('conn_mat.s',subs{k});
	##     eval([fieldname,'=fdt_network_matrix_norm'])
	## end
	## save('all_subj5000.mat',conn_mat)

# Results can be used in structural graph theory analyses
# See DTIConn.m
#=====================7a====================#


#=====================7b====================#
#======FOR USE WITH PROBTRACKX: N_SEED======#
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