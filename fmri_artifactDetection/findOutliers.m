%Script will output # of outliers generated by ART toolbox

%Load subject IDs
%should be in format: subs={'S03000','S03002','S03005','S03006'};
wave=load(subs.mat);

%Specify study directory
studyDir = '';

art_files={};
outliers={};

for j = 1:length(wave)
    subID=wave{j};
    try
        cd(strcat(studyDir,(subID),'/Faces/preprocess/'));
        art_files{j}=dir('art_regression_outliers_af*.mat');
        load(art_files{1,j}.name);
        outliers{j}=size(R,2);
    catch ME
        disp(['CHECK SUBJECT ',num2str(subID)])
        continue
    end
end

outliers={wave;outliers};