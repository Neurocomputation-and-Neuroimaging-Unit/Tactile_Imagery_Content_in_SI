function D6_glm_2ndLevel_OneSampleTTest(src_dir, SJin, inputfile, explicit_mask)

% This function performs a second level statistical analysis by performing
% one-sample t-tests over input files

[fdir fname fext]=fileparts(inputfile);
% target directory that will contain the created job.mat file
tgt_dir      = [src_dir filesep '2nd_level_' fdir];

% Create tgt_dir
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

% obtain the single subject image filenames
fNames = [];
for sj = 1:numel(SJin)
    temp_dir = fullfile(src_dir, SJin{sj}, fdir);
    fNames   = [fNames; cellstr([temp_dir filesep fname fext ',1']) ];
end

% create one sample t-test GLM
% -------------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir                       = {tgt_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = cellstr(fNames);
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = explicit_mask;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

%%  Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% T-Contrast Specification
% one-sample t-test contrast specification

% Allocate SPM.mat file
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(tgt_dir, 'SPM.mat')};
% Cycle over contrast specifications
% Allocate t-contrast structure
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'Positive Contrast';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% Delete existing contrasts (1=yes)
matlabbatch{3}.spm.stats.con.delete = 0;

%% Run the job ['Con' num2str(con) '_' cnames{con}]
fprintf(['Computing 2nd Level Contrast ' fdir '\n'])
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch
