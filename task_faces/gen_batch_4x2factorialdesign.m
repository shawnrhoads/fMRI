% Load subIDs 
% Should be in format: subID={'50016','50018','50022','50026'};
subID=load('subs.mat');

% Set studyDir
studyDir = '';

matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(studyDir)};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subjects';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'AttnVSIntrospect'; % set contrast name
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Face';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

for i=1:117
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {[studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0001.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0002.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0003.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0004.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0005.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0006.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0007.img,1']
        [studyDir,subID{i},'/Faces/1stlevel_contrasts/encoding/beta_0008.img,1']};
        
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).conds = [i 1 1
                                                                                      i 1 2
                                                                                      i 1 3
                                                                                      i 1 4
                                                                                      i 2 1
                                                                                      i 2 2
                                                                                      i 2 3
                                                                                      i 2 4];
end
 
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 3;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{3}.inter.fnums = [2
                                                                                  3];
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{4}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
