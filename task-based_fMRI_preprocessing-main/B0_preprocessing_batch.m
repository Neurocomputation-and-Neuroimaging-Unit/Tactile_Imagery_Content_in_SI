% function B0_preprocessing_batch_HIVR_tn

% please send questions to Till Nierhaus (till.nierhaus@fu-berlin.de) or Timo T. Schmidt (timo.t.schmidt@fu-berlin.de).

%### step A, structure data before running the preprocessing
%### --> convert DICOM images to 4D NIFTI image
%### --> subject folders (e.g. 'SJ002') including runs (e.g. 'run01') and anatomy (e.g. 'T1')
%############################################################################################

% required toolboxes:
% SPM12 ( http://www.fil.ion.ucl.ac.uk/spm/ )
% Rest_plus (1.8) ( http://restfmri.net/forum/index.php?q=rest )
% Step 7  Scrubbing: BRAMILA toolbox ( https://users.aalto.fi/~eglerean/bramila.html )
% Step 8  CompCorr: dPABI toolbox ( http://rfmri.org/dpabi )
% Step 11 Detrending: we use the LSGM-function ( as described in 10.1016/j.neuroimage.2003.12.042 )

clc
clear all
close all
addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data_processing_sw')

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
%SJs          = { 'SJ001' 'SJ002' 'SJ003' 'SJ004' 'SJ005' 'SJ006' 'SJ007' 'SJ008' 'SJ009' 'SJ010',...
%                 'SJ011' 'SJ012' 'SJ013' 'SJ014' 'SJ015' 'SJ016' 'SJ017' 'SJ018' 'SJ019' 'SJ020',...
%                 'SJ021' }; %  'SJ001' 'SJ004'

%welche SJ?
 % 1 2 4 5 6 7 8 % 10 12 13 14 16 17 18 % 19 20 21 22 23 24 25 % 26 27 28 29 30 31 32
sbjs = [26 27 28 29 30 31 32];
% nicht [3, 9, 11, 15]) %unvollständige datensätze

%sessions identifiers
runs = {'run01';'run02';'run03';'run04';'run05';'run06'};

%anatomy identifier
ana='T1';

%anatomical masks
wm_mask=['C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\wm_mask_eroded.nii']; %white matter mask file
csf_mask=['C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\csf_mask_eroded.nii']; %csf mask file
full_brain_mask=['C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\full_brain_mask.nii']; %full brain mask file

% selection of analysis steps (1-12) to be performed
analysis_switch = [11]; %1 4 3 5 7 8 9
start_prefix='wragf4d'; %'gf4d'; %ohne artrepair 


%# step 1  Segmentation
%# ------ Create nuisance masks on your own or take the provided ones
%# ------ NOTE: spike removal (e.g. "artrepair") should be performed as first step
%# step 2 --> remove first x scans                       --> prefix: x(number of cut volumes)
x=0;
%# step 3 --> slice time correction                      --> prefix: a
%  for interleaved slice order: do slice time correction, then realignment
%  otherwise do first realignment, then slice time correction (in analysis_switch 4 before 3)
n_slices = 37; % number of slices
slice_order=[1:n_slices];
refslice=19; % reference slice
TR=2; % repetition time in sec.
%# step 4  Realignment                                --> prefix: r
%# step 5  Coregister (estimate) mean-epi 2 anatomy
%# step 6  Coregister (estimate & resclice) mean-epi 2 anatomy --> prefix c
%# step 7  Normalization                              --> prefix: w
vox_size=[2 2 2]; % voxel size in mm
%# step 8  Scrubbing: calculate, interpolate outliers --> prefix: m(scrub_thresh)
scrub_thresh=0.4; % threshhold FD for scrubbing
%# step 9 Calculate WM and CSF Nuisance Signal
numComp = 5; % number of principle components
%# step 10 Smoothing                                   --> prefix: s
kernel_size=[3 3 3]; %FWHM kernel size
%# step 11 Detrending                                 --> prefix: d

%#####################################################
%#################### INPUT end ######################
%#####################################################

%%
currPrefix=start_prefix;

for n=analysis_switch
    
    switch n
        
        case 1 %% Segmentation
            warning off
            for sj = sbjs %sbjs %1:numel(SJs)
                if exist([src_dir filesep SJs{sj} filesep ana])==7
                    display(['Step 1, segmentation: ' SJs{sj}])
                    struct_dir = fullfile(src_dir, SJs{sj}, ana);
                    B1_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                else
                    display('###########################################################')
                    display(['############### ' SJs{sj} ', ' ana ' does not exsist ###########'])
                end
            end
            
        case 2 %% Delete first X scans
            if x>0
                for sj = sbjs %sbjs %1:numel(SJs)
                    for r = 1:numel(runs)
                        if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                            display(['Step 2, delete first ' num2str(x) ' volumes: ' SJs{sj} ', ' runs{r}])
                            run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                            B2_delete_scans(run_dir, ['^' currPrefix '.*\.nii'],x);
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                        end
                    end
                end
                currPrefix=['x' num2str(x) currPrefix];
            end
            
        case 3 %% Slice time correction
            for sj = sbjs %sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 3, slice time correction: ' SJs{sj} ', ' runs{r}])
                        run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        B3_slice_time_correction(SJs{sj},runs{r}, run_dir, ['^' currPrefix '.*\.nii'],n_slices,slice_order,refslice,TR);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['a' currPrefix];
            
        case 4 %% Realignment
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 4, realignment: ' SJs{sj} ', ' runs{r}])
                        run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        B4_realignment_run(run_dir, SJs{sj}, ['^' currPrefix '.*\.nii']);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['r' currPrefix];
            
        case 5 %% Coregister (estimate) mean-epi 2 anatomy
            warning off
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 5, coregistration: ' SJs{sj} ', ' runs{r}])
                        func_dir        = fullfile(src_dir, SJs{sj}, runs{r});
                        struct_dir      = fullfile(src_dir, SJs{sj}, ana);
                        B5_coregister_est(currPrefix, func_dir, struct_dir, SJs{sj}, '^s.*\.nii');
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            
        case 6 %% Coregister (estimate & reslice) mean-epi 2 anatomy
            warning off
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 6, coregistration: ' SJs{sj} ', ' runs{r}])
                        func_dir        = fullfile(src_dir, SJs{sj}, runs{r});
                        struct_dir      = fullfile(src_dir, SJs{sj}, ana);
                        B6_coregister_est_re(currPrefix, func_dir, struct_dir, SJs{sj}, '^s.*\.nii');
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            
        case 7 %% Normalization
            for sj = sbjs %1:numel(SJs)
                struct_dir = fullfile(src_dir, SJs{sj}, ana);
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 7, normalization: ' SJs{sj} ', ' runs{r}])
                        data_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        B7_normalization_run(data_dir, struct_dir, SJs{sj}, ['^' currPrefix '.*\.nii'], vox_size);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['w' currPrefix];
            
        case 8 %% Scrubbing: calculate outliers
            scrub_prefix=['m' num2str(scrub_thresh)];
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 8, scrubbing: ' SJs{sj} ', ' runs{r}])
                        data_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        %estimate and save motion statistics
                        n=1;
                        f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) '.*\.txt']);
                        while isempty(f)
                            n=n+1;
                            f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) '.*\.txt']);
                        end
                        cfg.motionparam=[data_dir filesep f];
                        cfg.prepro_suite = 'spm';
                        [fwd,rms]=bramila_framewiseDisplacement(cfg);
                        outliers=fwd>scrub_thresh;
                        percent_out=(sum(outliers)/length(outliers))*100;
                        disp(['outliers for ' num2str(SJs{sj}) ', ' runs{r} ': ' num2str(percent_out) '%']);
                        save([data_dir filesep scrub_prefix currPrefix '_' SJs{sj} '_' runs{r} '_FWDstat.mat'],'fwd','rms','outliers','percent_out','scrub_thresh','cfg')
                        %srub outliers by replacing them with average of nearest neighbors
                        B8_scrub_data(data_dir, ['^' currPrefix '.*\.nii'], outliers,  scrub_prefix);
                        all_percent_out(sj,r)=percent_out;
                        all_rp{sj,r}=load(cfg.motionparam);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=[scrub_prefix currPrefix];
            save([src_dir filesep 'all_MOTIONstat_' currPrefix '.mat'],'SJs','runs','scrub_thresh','all_percent_out','all_rp')
            
        case 9 %% CompCorr
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 9, CompCorr: ' SJs{sj} ', ' runs{r}])
                        data_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        B9_compcorr_run(data_dir, SJs{sj}, ['^' currPrefix '.*\.nii'], numComp, wm_mask, csf_mask);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            
        case 10 %% Smoothing
            for sj = sbjs %1:numel(SJs)
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        display(['Step 10, smoothing: ' SJs{sj} ', ' runs{r}])
                        run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        B10_smoothing_run(run_dir, SJs{sj}, ['^' currPrefix '.*\.nii'],kernel_size);
                        display([SJs{sj} ', ' runs{r} ' is done'])
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['s' num2str(unique(kernel_size)) currPrefix];
            
        case 11 %% Detrending
            for sj = sbjs %1:numel(SJs)
                display(['Step 11, Detrending: ' SJs{sj}])
                for r = 1:numel(runs)
                    if exist([src_dir filesep SJs{sj} filesep runs{r}])==7
                        run_dir = fullfile(src_dir, SJs{sj}, runs{r});
                        f = spm_select('List',run_dir, ['^' currPrefix '.*\.nii']);
                        V=spm_vol([run_dir filesep f(1,:)]);
                        files={};
                        for i=1:(size(V,1))
                            files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
                        end
                        fileset{r}=char(files);
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{r} ' does not exsist ###########'])
                    end
                end
                if exist('fileset')
                    B11_detrending_lmgs(fileset);
                    clear fileset files
                end
                display([SJs{sj} ', is done'])
            end
        
        otherwise
            display('######################################################################################')
            display(['############################## Case ' num2str(n) ' does not exsist ##############################'])
            display('######################################################################################')
    end
end
%%
% adjust the following to delete files with certain prefixes from the subject
% folders
if ismember(99,analysis_switch)
    for sj = sbjs %sbjs %1:numel(SJs)
        for runnr = 1:numel(runs)
            data_dir = fullfile(src_dir, SJs{sj}, ['run0', num2str(runnr)]);
            cd(data_dir)
            delete('r*.mat')
            delete('f*.mat')
            delete('wraf4d*.mat')
            delete('ds8wraf4d*.mat')
            delete('ds8m0.4wraf4d*.mat')
            delete('m*.mat')
            delete('a*.mat')
%             delete('o*.mat')
%             delete('d*.mat')
%             delete('s*.mat')
%             delete('w*.mat')
        end
    end
end