# misc_brain contains scripts for fMRI/DTI preprocessing and analysis

<b>./fMRI_artifactDetection</b> contains two scripts for recursively generating configuration settings per subject in study and running ART toolbox

<b>./probTractography</b> contains script to recursively pre-process diffusion tensor imaging data and run probabilistic tractography analysis. (masterFDT_parallel.sh runs subjects in parallel).
Script can use two methods:
    1) To generate a NxN connectivity matrix;
    2) To run a seed2target analysis

<b>./task_faces</b> contains 2nd level analysis batch script for 4x2 factorial design

<b>./task_whichFear</b> to be updated soon