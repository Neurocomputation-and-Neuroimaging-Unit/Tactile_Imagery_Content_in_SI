% function D0_decoding_batch_HIVR

% please send questions to Till Nierhaus (till.nierhaus@fu-berlin.de) or Timo T. Schmidt ().

%### step A, structure data before running the preprocessing
%### --> convert DICOM images to 4D NIFTI image
%### --> subject folders (e.g. 'SJ002') including runs (e.g. 'run01') and anatomy (e.g. 'T1')
%############################################################################################

% required toolboxes:
% the decoding toolbox TDT
%  https://doi.org/10.3389/fninf.2014.00088

clc
clear all
close all
addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data_processing_sw')
addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\spm12')
addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\decoding_toolbox')

%#####################################################
%#################### INPUT ##########################
%#####################################################

%SPM-path
SPM_path  = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\spm12';

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data';

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end
% SJs          = { 'SJ001' 'SJ002' 'SJ003' 'SJ004' 'SJ005' 'SJ006' 'SJ007' 'SJ008' 'SJ009' 'SJ010',...
%                  'SJ011' 'SJ012' 'SJ013' 'SJ014' 'SJ015' 'SJ016' 'SJ017' 'SJ018' 'SJ019' 'SJ020',...
%                  'SJ021' }; %  'SJ001' 'SJ004'

% sbjs = [1 2 4 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32]; 
sbjs = [1 2 6 7 8 10 12 13 14 16 18 19 20 21 22 24 26 27 29 30 32]; 
% 1 2 4 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
% ohne 25 weil ein run wseniger und das ist nicht so smart

%sessions identifiers
runs = {'run01';'run02';'run03';'run04';'run05';'run06'};

%anatomy identifier
ana='T1';

% selection of analysis steps (1-5) to be performed
analysis_switch = [4]; %2 1 3 4 5
start_prefix='wragf4d'; %'f4d'; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1  --> slice time correction                      --> prefix: a
%  for interleaved slice order: do slice time correction, then realignment
%  otherwise do first realignment, then slice time correction (in analysis_switch 2 before 1)
n_slices = 37; % number of slices
slice_order=[1:n_slices];
refslice=19; % reference slice
TR=2; % repetition time in sec.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  Realignment for all runs!!!                   --> prefix: r2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 3: 1st level glm for NOT-normalized data
log_folder = 'logs'; % folder with logfiles inside subject directory
logfile_prefix = 'all_onsets_goodImag_'; % file prefix of .mat file for each subject;
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec

%sessions identifiers
condnames = {'StimPress','StimFlutt','StimVibro','ImagPress','ImagFlutt','ImagVibro','Null_1','Null_2','preCue','Response','badImag'};  % condition names
duration = 3; % epoch duration; for single events set 0
tr   = 2.0;     % TR
fmri_t    = 37; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0   = 19; % Microtime onset; If you have performed slice-timing correction, you must change this parameter to match the reference slice specified there; otherwise 1 for first, fmri_t for last slice, or fmri_t/2 as compromise
hpf      = 128; % High-pass filter cut-off; default 128 sec
% include multiple regressors (1=yes)
hm=0;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=0;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 4: Decoding
%'SET THE DIRECTORY WHERE YOUR SPM.MAT AND ALL RELATED BETAS ARE, individual 1st level glm directory';
betaDir4 = '1st_level_D0_preCue_goodImag_wragf4d';
%'SPECIFY WHERE THE RESULTING ACCURACY MAPS SHOULD BE STORED';
outputDir4 = 'rev_wragf4d_1C_preCue_goodImag_maskCode425';
% Set the label name pairs to the regressor names which you want to use for decoding, e.g. 'higher' and 'lower'
% if you don't remember the names --> run display_regressor_names(beta_dir)
labelnames = {'StimPress','StimFlutt', 'StimVibro'};
theseConds = '_Stims';
% {'StimPress','StimFlutt', 'StimVibro', 'ImagPress','ImagFlutt','ImagVibro'};
% {'ImagPress','ImagFlutt','ImagVibro'};
% {'StimPress','StimFlutt', 'StimVibro'};

%     {'StimPress','StimFlutt'; ...
%     'StimPress','StimVibro'; ...
%     'StimFlutt','StimVibro'; ...
%     'ImagPress','ImagFlutt'; ...
%     'ImagPress','ImagVibro'; ...
%     'ImagFlutt','ImagVibro'}; 
%     'ImagPress','Att'; ...
%     'ImagFlutt','Att'; ...
%     'ImagVibro','Att'};
%     'StimPress','Att'; ...
%     'StimFlutt','Att'; ...
%     'StimVibro','Att'
%     'StimPress','ImagPress'; ...
%     'StimFlutt','ImagFlutt'; ...
%     'StimVibro','ImagVibro'};
excep_array = [];
% if strcmp(string(labelnames(end,2)), 'Att')
%     excep_array = size(labelnames,1)-2:size(labelnames,1);
% else
%     excep_array = [];
% end
% rel_files = {strcat('ImagPress_Att\', 's5wres_accuracy_minus_chance.nii'), strcat('ImagFlutt_Att\', 's5wres_accuracy_minus_chance.nii'), strcat('ImagVibro_Att\', 's5wres_accuracy_minus_chance.nii')};
% rel_files = {strcat('StimPress_StimFlutt\', 's5wres_accuracy_minus_chance.nii'), ...
%              strcat('StimPress_StimVibro\', 's5wres_accuracy_minus_chance.nii'), ...
%              strcat('StimFlutt_StimVibro\', 's5wres_accuracy_minus_chance.nii')}%, ...
%              strcat('ImagPress_ImagFlutt\', 's5wres_accuracy_minus_chance.nii'), ...
%              strcat('ImagPress_ImagVibro\', 's5wres_accuracy_minus_chance.nii'), ...
%              strcat('ImagFlutt_ImagVibro\', 's5wres_accuracy_minus_chance.nii')};

norm4 = 0; %normalise accuracy maps (1=yes)
vox_size=[2 2 2]; % voxel size in mm
smoo4 = 0; %smooth with kernel in mm
ostt4 = 0; % 2nd-level one-sample ttest for all label-pairs (1=yes)
exclSJs4 = [3 9 11 15 25]; % subjects to be excluded from 2nd-level analysis
explicit_mask4={'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\brain_mask.nii'}; % specify explicit mask file for 2nd-level analysis, else {''}

cd('C:\Users\saraw\Desktop\newest_masks\GoodImag_model\controlClustersConjUncorr_k200');
roi_names = dir(['*.nii']);
% roi_names = 'C:\Users\saraw\Desktop\Masks\MASKS2\thatmask\rCONJ_att_CUT_allOfSomato.nii';

% SET OTHER DECODING CFG PARAMETERS in D2_Decoding_batch.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 5: Decoding, cross classification
%'SET THE DIRECTORY WHERE YOUR SPM.MAT AND ALL RELATED BETAS ARE, individual 1st level glm directory';
betaDir5 = '1st_level_D0_preCue_wragf4d';
%'SPECIFY WHERE THE RESULTING ACCURACY MAPS SHOULD BE STORED';
outputDir5 = 'rev_XClass_wragf4d_1C_preCue_goodImag_maskCode425';
% Set the label name pairess to the regressor names which you want to use for decoding, e.g. 'higher' and 'lower'
% if you don't remember the names --> run display_regressor_names(betaDir)
% labels of training and test data must match for proper cross classification! 
labelnames_train = {'StimPress','StimFlutt', 'StimVibro'};

% {'StimPress','StimFlutt'; ...
%     'StimPress','StimVibro'; ...
%     'StimFlutt','StimVibro';...
%     'ImagPress','ImagFlutt'; ...
%     'ImagPress','ImagVibro'; ...
%     'ImagFlutt','ImagVibro'};

labelnames_test  = {'ImagPress','ImagFlutt', 'ImagVibro'};

% {'ImagPress','ImagFlutt'; ...
%     'ImagPress','ImagVibro';...
%     'ImagFlutt','ImagVibro'; ...
%     'StimPress','StimFlutt'; ...
%     'StimPress','StimVibro'; ...
%     'StimFlutt','StimVibro'};

% roi_names = {['PSC_1_cut_interBA_cut_ROI_left_SI_MNI_cut_STIM.nii'], ...
%              ['PSC_1_cut_interBA_cut_ROI_right_SI_MNI_cut_STIM.nii'], ...
%              ['PSC_2_cut_interBA_cut_ROI_left_SI_MNI_cut_STIM.nii'], ...
%              ['PSC_2_cut_interBA_cut_ROI_right_SI_MNI_cut_STIM.nii'], ...
%              ['PSC_3b_cut_interBA_cut_ROI_left_SI_MNI_cut_STIM.nii'], ...
%              ['PSC_3b_cut_interBA_cut_ROI_right_SI_MNI_cut_STIM.nii'], ...
%              ['ROI_SIILeftMask_MNI_cut_NEW_StimMask001.nii'], ...
%              ['ROI_SIIRightMask_MNI_cut_NEW_StimMask001.nii']};
%          
% dec2_files = {strcat('TRAIN_StimPress_StimFlutt_TEST_ImagPress_ImagFlutt\', 'res_accuracy_minus_chance'), ...
%              strcat('TRAIN_StimPress_StimVibro_TEST_ImagPress_ImagVibro\', 'res_accuracy_minus_chance'), ...
%              strcat('TRAIN_StimFlutt_StimVibro_TEST_ImagFlutt_ImagVibro\', 'res_accuracy_minus_chance'), ...
%              strcat('TRAIN_ImagPress_ImagFlutt_TEST_StimPress_StimFlutt\', 'res_accuracy_minus_chance'), ...
%              strcat('TRAIN_ImagPress_ImagVibro_TEST_StimPress_StimVibro\', 'res_accuracy_minus_chance'), ...
%              strcat('TRAIN_ImagFlutt_ImagVibro_TEST_StimFlutt_StimVibro\', 'res_accuracy_minus_chance')};

norm5 = 0; %normalise accuracy maps (1=yes)
vox_size=[2 2 2]; % voxel size in mm
smoo5 = 0; %smooth with kernel in mm
ostt5 = 0; % 2nd-level one-sample ttest for all label-pairs (1=yes)
exclSJs5 = [3 9 11 15 25]; % subjects to be excluded from 2nd-level analysis
explicit_mask5={'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\brain_mask.nii'}; % specify explicit mask file for 2nd-level analysis, else {''}

% SET OTHER DECODING CFG PARAMETERS in D3_Decoding_batch_XClass.m


%#####################################################
%#################### INPUT end ######################
%#####################################################

%%
currPrefix=start_prefix;

for n=analysis_switch
    
    switch n
        
        case 1 %% Slice time correction
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 1, slice time correction: ' SJs{sj} ', ' runs{r}])
                        run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        D1_slice_time_correction(SJs{sj},runs{r}, run_dir, ['^' currPrefix '.*\.nii'],n_slices,slice_order,refslice,TR);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['a' currPrefix];
            
        case 2 %% Realignment all runs
            for sj = sbjs %1:numel(SJs)
                display(['Step 2, realignment all runs: ' SJs{sj} ])
                sj_dir = fullfile(src_dir, SJs{sj});
                D2_Realignment_all_runs(sj_dir, runs, ['^' currPrefix '.*\.nii']);
            end
            currPrefix=['r2' currPrefix];
            
        %% ab hier
            
        case 3 %% step 3: 1st level glm, cycle over subjects
            
            beta_dir = ['1st_level_D0_preCue_goodImag_' currPrefix]; % folder that will contain the created job.mat file and SPM file
            
            for sj = sbjs %1:numel(SJs)
                display(['Step 3, 1st level glm: ' SJs{sj} ])
                subj_dir = fullfile(src_dir, SJs{sj});
                try
                    load([subj_dir filesep log_folder filesep logfile_prefix SJs{sj} '.mat'],'onsets');
                    D3_glm_1stLevel(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc);
                catch
                    display('###################################################################')
                    display(['################## ' SJs{sj} ', ERROR GLM ###################'])
                    display('###################################################################')
                end
            end
            
        %% und hier
            
        case 4 %% step 4: Decoding

%             parameterC = [];
            
            for sj = sbjs %1:numel(SJs)
%                 try
                    for Lp=1:size(labelnames,1)
                        
                        Lp

                        s = find(sbjs == sj);

                        display(['Step 4, Decoding: ' SJs{sj} ', Conditions: ' cell2mat(labelnames(Lp,1)) ', ' cell2mat(labelnames(Lp,2))])
                        beta_path    = fullfile(src_dir, SJs{sj}, betaDir4);
%                         output_path  = fullfile(src_dir, SJs{sj}, [outputDir4 '_' cell2mat(labelnames(Lp,1)) '_' cell2mat(labelnames(Lp,2)) '_' cell2mat(labelnames(Lp,3))]); %% vlt mal drüber nachdenken, noch mehr unterordner zu machen
                                                output_path  = fullfile(src_dir, SJs{sj}, [outputDir4 theseConds]); %% vlt mal drüber nachdenken, noch mehr unterordner zu machen

                        D4_Decoding_batch(beta_path, output_path, labelnames(Lp,:)) %% hier is das wichtige
%                         
                        Images2 = {};
                        
                        if norm4
                            Images1{Lp}=fullfile(output_path,'res_accuracy_minus_chance.nii');
                            if smoo4
                                Images2{Lp}=fullfile(output_path,'wres_accuracy_minus_chance.nii');
                                if ~ismember(Lp, excep_array)
                                    I4ostt{Lp}=fullfile([outputDir4 '_Imags'],['s' num2str(smoo4) 'wres_accuracy_minus_chance.nii']);
                                end
                            else
                                if ~ismember(Lp, excep_array)
                                    I4ostt{Lp}=fullfile([outputDir4 '_' cell2mat(labelnames(Lp,1)) '_' cell2mat(labelnames(Lp,2))],'wres_accuracy_minus_chance.nii');
                                end
                            end
                        elseif smoo4
%                             for roi = 1:length(roi_names)
%                                 Images2=[Images2; fullfile(output_path, strcat('res_accuracy_minus_chance_', roi_names(roi)))];
%                                 I4ostt=[I4ostt; fullfile(strcat(outputDir4, '_', cell2mat(labelnames(Lp,1)), '_', cell2mat(labelnames(Lp,2)), 's', num2str(smoo4), 'wres_accuracy_minus_chance_', roi_names(roi)))];
%                             end
                            Images2{Lp}=fullfile(output_path,'res_accuracy_minus_chance.nii');

                        end
                        
                        close all
                    end

                    if norm4
                            func_dir    = fullfile(src_dir, SJs{sj}, runs{1});
                            struct_dir  = fullfile(src_dir, SJs{sj}, ana);
                            D4a_coregister_est(func_dir, struct_dir, Images1)
                            fprintf(['Normalizing accuracy maps ', SJs{sj}, '\n']);
                            D4b_normalization(Images1, struct_dir, vox_size)
                    end
                    if smoo4
                        fprintf(['Smoothing accuracy maps ', SJs{sj}, '\n']);
                        D4c_smoothing(Images2, smoo4)
                    end

                    if length(excep_array) > 0 %%%%%%
                        
                        matlabbatch{1}.spm.util.imcalc.input = {strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_', rel_files{1})
                                                                strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_', rel_files{2})
                                                                strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_', rel_files{3})};

%                                 matlabbatch{1}.spm.util.imcalc.output = 'res_accuracy_minus_chance.nii';

                        matlabbatch{1}.spm.util.imcalc.output = 's5wres_accuracy_minus_chance.nii';
                        
                        if ~exist(strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_meanStim'), 'dir')
                            mkdir(strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_meanStim'));
                        end
                        
                        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(src_dir, filesep, SJs{sj}, filesep, outputDir4, '_meanStim')};

                        matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 + i3) / 3';

                        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

                        I4ostt{Lp}=fullfile([outputDir4 '_meanStim'],'s5wres_accuracy_minus_chance.nii');
                        I4ostt(excep_array(1:2)) = []

                        spm('defaults', 'FMRI');
                        spm_jobman('run', matlabbatch);
                    
                    end
                    clear Images*
%                 catch
%                     display('###################################################################')
%                     display(['################### ' SJs{sj} ', ERROR DECODING ###################'])
%                     display('###################################################################')
%                 end
            end
            
            if ostt4
                SJin=SJs;
                SJin(exclSJs4)=[];
                for t=1:size(I4ostt,2)
                    D6_glm_2ndLevel_OneSampleTTest(src_dir, SJin, I4ostt{t}, explicit_mask4)
                end
            end
            clear I4ostt SJin
            
        %% und später hier
            
        case 5 %% step 5: Decoding, cross classification
            if size(labelnames_train)==size(labelnames_test)
                
                for sj = sbjs %1:numel(SJs)
    %                 try
                        for Lp=1:size(labelnames_train,1)
                            display(['Step 5, Decoding XClass: ' SJs{sj} ]) 
                            beta_path    = fullfile(src_dir, SJs{sj}, betaDir5);
%                             output_path  = fullfile(src_dir, SJs{sj}, [outputDir5 '_TRAIN_' cell2mat(labelnames_train(Lp,1)) '_' cell2mat(labelnames_train(Lp,2)) '_TEST_' cell2mat(labelnames_test(Lp,1)) '_' cell2mat(labelnames_test(Lp,2))]);
                            output_path  = fullfile(src_dir, SJs{sj}, [outputDir5 '_TRAIN_Stims_TEST_Imags']);
                            D5_Decoding_batch_XClass(beta_path, output_path, labelnames_train(Lp,:), labelnames_test(Lp,:))% hier is das andere wichtige

                            Images2 = {};
                            I5ostt = {};

                            if norm5
                                Images1{Lp}=fullfile(output_path,'res_accuracy_minus_chance.nii');
                                if smoo5
                                    Images2{Lp}=fullfile(output_path,'wres_accuracy_minus_chance.nii');
                                    I5ostt{Lp}=fullfile([outputDir5 '_TRAIN_Stims_TEST_Imags'],['s' num2str(smoo5) 'wres_accuracy_minus_chance.nii']);
                                else
                                    I5ostt{Lp}=fullfile([outputDir5 '_TRAIN_Stims_TEST_Imags'],'wres_accuracy_minus_chance.nii');
                                end
                            elseif smoo5
                                for roi = 1:length(roi_names)
                                    Images1=[Images1; fullfile(output_path, strcat('res_accuracy_minus_chance_', roi_names(roi)))];
                                    I5ostt=[I4ostt; fullfile(strcat(outputDir5, ['_TRAIN_' cell2mat(labelnames_train(Lp,1)) '_' cell2mat(labelnames_train(Lp,2)) '_TEST_' cell2mat(labelnames_test(Lp,1)) '_' cell2mat(labelnames_test(Lp,2))], 's', num2str(smoo5), 'wres_accuracy_minus_chance_', roi_names(roi)))];
                                end
    %                             Images2{Lp}=fullfile(output_path,'res_accuracy_minus_chance.nii');
                            end
                        end
                        close all

                        if norm5
                            func_dir    = fullfile(src_dir, SJs{sj}, runs{1});
                            struct_dir  = fullfile(src_dir, SJs{sj}, ana);
                            D4a_coregister_est(func_dir, struct_dir, Images1)
                            fprintf(['Normalizing accuracy maps ', SJs{sj}, '\n']);
                            D4b_normalization(Images1, struct_dir, vox_size)
                        end
                        if smoo5
                            fprintf(['Smoothing accuracy maps ', SJs{sj}, '\n']);
                            D4c_smoothing(Images2, smoo5)
                        end

                        clear Images*

%                         end
    %                 catch
    %                     display('###################################################################')
    %                     display(['########### ' SJs{sj} ', ERROR DECODING X Classification ###########'])
    %                     display('###################################################################')
    %                 end
                end

                if ostt5
                    SJin=SJs;
                    SJin(exclSJs5)=[];
                    for t=1:size(I5ostt,2)
                        D6_glm_2ndLevel_OneSampleTTest(src_dir, SJin, I5ostt{t}, explicit_mask5)
                    end
                end
                clear I5ostt SJin

            else
                display('###################################################################')
                display('######################## Input Error ##############################')
                display('############## Training and Test Labels must match! ###############')
                display('###################################################################')
            end
    end
end
%%
% adjust the following to delete files with certain prefixes from the subject
% folders
% if ismember(99,analysis_switch)
%     for sj = 1:numel(SJs)
%         for runnr = 1:numel(runs)
%             data_dir = fullfile(src_dir, SJs{sj}, ['run0', num2str(runs(runnr))]);
%             cd(data_dir)
%             delete('hp*.mat')
%             delete('o*.mat')
%             delete('d*.mat')
%             delete('s*.mat')
%             delete('w*.mat')
%         end
%     end
% end