%This script takes an array of subject numbers, stores all image names and
%locations as one string, each separated by one white space. 
%Then, uses this line to create an individual .cfg configuration file for
%that specific subject and saves it to '/usr/local/art'.
%
%The last set of lines can be used to copy and paste into the command line
%in case there are different files names for images and/or motion parameter
%files.
%
%The script also assumes that each image and motion parameter file share
%the same pattern for naming.. It can be edited to account for differences
%in naming.

addpath /usr/local/art
addpath /usr/local/spm8

%Load subject IDs
%should be in format: subs={'S03000','S03002','S03005','S03006'};
subs=load(subs.mat);

%Specify study directory
studyDir = '';

%Specify location of slice-timing corrected functional niftys 
STdir = '';

%Specify total number of volumes
totVol = ;

cfg_lines={};

%Start generating a new cfg file for each subject in '/usr/local/ART'
for j = 1:length(subs)
    subID=subs{j};
    subID_corr=subs_corrected{j};
    for i = 1:totVol
        cfg_lines{i}=fullfile(studydir,num2str(subID),STdir,'af',num2str(subID),'*',num2str(i),'.nii'];
    end
    
    cfg_1line = strjoin(cfg_lines,' ');

%     Un-comment if you would like to save the .mat for each subject in current folder.
%     filename = ['ARTcfg/cfg',num2str(subID),'.mat'];
%     save(filename,'cfg_1line')

    %Write new cfg script for ART
    fid = fopen(['/usr/local/art/',num2str(subID),'_config.cfg'],'wt');
    fprintf(fid, 'sessions: 1\n');
    fprintf(fid, 'global_mean: 1\n');
    fprintf(fid, 'global_threshold: 3.0\n');
    fprintf(fid, 'motion_threshold: 0.5\n');
    fprintf(fid, 'motion_file_type: 0\n');
    fprintf(fid, 'motion_fname_from_image_fname: 0\n');
    fprintf(fid, 'end\n');
    fprintf(fid, ['session 1 image ',cfg_1line]);
    fprintf(fid, '\n');
    
    fprintf(fid, ['session 1 motion ',studyDir,num2str(subID),STdir,'rp_af',num2str(subID),'-0007-00001-000001-01.txt\n']);
    fprintf(fid, 'end\n');
    fclose(fid);
end


%Run ART
count=0;
for k=1:length(subs)
    count=count+1
    percent=(count/length(subs))*100;
    disp([num2str(percent),'%'])
    subID=subs{k};
    subID_corr=subs_corr{k};
    try
        subDir=strcat(studyDir,num2str(subID),STdir);
        cd(subDir);
        art('sess_file',[subID_corr])
    catch ART_err
        disp(['CHECK SUBJECT ',num2str(subID)])
        continue
    end
end

%to replace text using bash: 'replace "d: 2.0" "d: 0.5" -- *_config.cfg'
