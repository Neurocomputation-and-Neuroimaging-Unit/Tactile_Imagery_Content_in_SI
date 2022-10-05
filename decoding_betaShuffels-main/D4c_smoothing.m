function D4c_smoothing(Images, ks)

warning off

kernel_size=[ks ks ks];

matlabbatch{1}.spm.spatial.smooth.data = Images';
matlabbatch{1}.spm.spatial.smooth.fwhm = kernel_size;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(ks)];

inputs = cell(0, 1);

% save the job variable to disc
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch, '', inputs{:});
% spm_jobman('serial', matlabbatch, inputs{:});
clear jobs

