%% create ROIs from Conjunction x ROIs_AnatomyToolbox

% ROIs AT

% threshold, unilat, dimensions, composite SII, cut

% volumeSize = squeeze(that.D);





clear anat_names anat
anat_path = 'C:\Users\saraw\Desktop\Masks\anatomy\';
cd(anat_path);
anat_names = dir('PSC*');
% roi_names([1 2], :) = [];
for an=1:length(anat_names)
    anat{an} = [anat_path anat_names(an).name];
end
tgt_dir = 'C:\Users\saraw\Desktop\Masks\TR40';
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

count_ind = [1 2 3; 2 1 3; 3 1 2];

for i = 1:size(anat_names, 1)

    clear thisANA

    thisANA = anat_names(i).name;
    length_name = length(thisANA);
    thisANA = thisANA(1:length_name-4);

    matlabbatch{1}.spm.util.imcalc.input = {strcat(anat{count_ind(i,1)}, ',1')
        strcat(anat{count_ind(i,2)}, ',1')
        strcat(anat{count_ind(i,3)}, ',1')};

    matlabbatch{1}.spm.util.imcalc.output = strcat(thisANA, '_TR40');

    matlabbatch{1}.spm.util.imcalc.outdir = {'C:\Users\saraw\Desktop\Masks\TR40\'};

    matlabbatch{1}.spm.util.imcalc.expression = '((i1 > 0.4) > i2) & ((i1 > 0.4) > i3)';
    %                 matlabbatch{1}.spm.util.imcalc.expression = 'i1 > 0.5';

    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

end


clear TR_names TRs
TR_path = 'C:\Users\saraw\Desktop\Masks\TR40\';
cd(TR_path);
TR_names = dir('PSC*');
for ma=1:length(TR_names)
    TRs{ma} = [TR_path TR_names(ma).name];
end

for tr = 1:size(TR_names, 1)

    clear thisTR

    thisTR = TR_names(tr).name;
    length_name = length(thisTR);
    thisTR = thisTR(1:length_name-4);

    for sides = [1 2]

        new_nii = load_untouch_nii(TRs{tr});

        if sides == 1

            side = 'left';

            new_nii.img([76:151],:,:) = 0;

        elseif sides == 2

            side = 'right';

            new_nii.img([1:75],:,:) = 0;

        end

        output = strcat(thisTR, '_', side);
        %                     info = niftiinfo(['C:\Users\saraw\Desktop\Masks\thesholded\' thisTR '.nii']);
        %                     info.Filename = ['C:\Users\saraw\Desktop\Masks\unilat\' output '.nii'];
        %                     info.ImageSize = [79, 95, 79];
        %                     info.PixelDimensions = [2 2 2];
        %                     info
        %                     new_nii = imresize3(new_nii, [79, 95, 79]);
        outdir = 'C:\Users\saraw\Desktop\Masks\unilat40\';
        if ~exist(outdir, 'dir')
            mkdir(outdir)
        end
        cd(outdir);
        save_untouch_nii(new_nii, output);

    end

end

clear uni_names rois
uni_path = 'C:\Users\saraw\Desktop\Masks\unilat50\realigned\SII\';
cd(uni_path);
uni_names = dir('*left.nii');
%             uni_names([1 2], :) = [];
for ma=1:length(uni_names)
    unis{ma} = [uni_path uni_names(ma).name];
end

%             for u = 1:size(uni_names, 1)
%
%                 clear thisROI
%
%                 thisROI = uni_names(u).name;
%                 length_name = length(thisROI);
%                 thisROI = thisROI(1:length_name-4);

matlabbatch{1}.spm.util.imcalc.input = {strcat(unis{1}, ',1')
    strcat(unis{2}, ',1')
    strcat(unis{3}, ',1')
    strcat(unis{4}, ',1')};

matlabbatch{1}.spm.util.imcalc.output = strcat('rSII_TR50_left');

matlabbatch{1}.spm.util.imcalc.outdir = {'C:\Users\saraw\Desktop\Masks\unilat50\realigned\'};

matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 + i3 + i4) > 0';

matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

%             end


clear roi_names rois
roi_path = 'C:\Users\saraw\Desktop\Masks\unilat\realigned\';
cd(roi_path);
roi_names = dir('r*');
%             roi_names([1 2], :) = [];
for ma=1:length(roi_names)
    rois{ma} = [roi_path roi_names(ma).name];
end

all_vs = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data\2nd_level_27_y_ds8wragf4d_FlexFact\Null_CONJ0001.nii';
outpath = 'C:\Users\saraw\Desktop\Masks\Null_CONJ\uncorr'

if ~exist(outpath, 'dir')
    mkdir(outpath)
end

for j = 1:size(rois, 2)

    clear thisROI

    thisROI = roi_names(j).name;
    length_name = length(thisROI);
    thisROI = thisROI(1:length_name-4);

    matlabbatch{1}.spm.util.imcalc.input = {
        strcat(all_vs, ',1')
        strcat(rois{j}, ',1')
        };

    matlabbatch{1}.spm.util.imcalc.output = strcat(thisROI, '_CUT_Null_CONJ0001');
    %                 matlabbatch{1}.spm.util.imcalc.output = strcat('SII_', j, '_CUT_all_vs_Att');


    matlabbatch{1}.spm.util.imcalc.outdir = {outpath};

    matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2) > 1.1';

    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

end


