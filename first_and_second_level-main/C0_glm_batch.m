% function glm_1stLevel_batch
%%
% This is a Batch-Script to compute GLM for all subjects
% The pre-processed data for each run is taken individually

clc
clear all
close all


%#####################################################
%#################### INPUT ##########################
%#####################################################

addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data_processing_sw')

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data';

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

sbjs = [1 2 4 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32]; % 1 2 4 5 6 7 8 % 10 12 13 14 16 17 18 % 19 20 21 22 23 24 % 26 27 28 29 30 31 32
%sbjs = [1 2 6 7 8 10 12 13 14 16 18 19 20 21 22 24 26 27 28 29 30 32]
%ohne 4 5 17 23 31 da sonst leere bedingungen
% SJs          = { 'SJ001' 'SJ002' 'SJ003' 'SJ004' 'SJ005' 'SJ006' 'SJ007' 'SJ008' 'SJ009' 'SJ010',...
%                  'SJ011' 'SJ012' 'SJ013' 'SJ014' 'SJ015' 'SJ016' 'SJ017' 'SJ018' 'SJ019' 'SJ020',...
%                  'SJ021' };  %  'SJ001' 'SJ004'

% selection of analysis steps to be performed
analysis_switch = [4];  % 1 2 3 4

%%sp
%# step 1: 1st level glm
log_folder = 'logs'; % folder with logfiles inside subject directory
logfile_prefix = 'all_onsets4_'; % file prefix of .mat file for each subject;
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec

%sessions identifiers
runs = {'run01';'run02';'run03';'run04';'run05';'run06'}; % 1 2 3 4 5 6
%specify filter for finding functional data
prefix_func = 'ds8wragf4d';

condnames = {'StimPress','StimFlutt','StimVibro','ImagPress','ImagFlutt','ImagVibro','Null_1','Null_2','Att_1','Att_2','Response'};  % condition names
duration = 3; % epoch duration; for single events set 0
outputfolder_1st = ['1st_level_27_z_' prefix_func]; % folder that will contain the created job.mat file and SPM file
tr   = 2.0;     % TR
fmri_t    = 37; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0   = 19; % Microtime onset; If you have performed slice-timing correction, you must change this parameter to match the reference slice specified there; otherwise 1 for first, fmri_t for last slice, or fmri_t/2 as compromise
hpf      = 128; % High-pass filter cut-off; default 128 sec

% include multiple regressors (1=yes)
hm=0;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=0;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
    % if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st
    
%%
%# step 2: contrasts for 1st level
analysisfolder = '1st_level_27_z_ds8wragf4d'; % folder that contains SPM file
% T-Contrast Specification
cnames = {'StimPress' , ...         %1
          'StimFlutt' , ...         %2
          'StimVibro' , ...         %3
          'ImagPress' , ...         %4
          'ImagFlutt' , ...         %5
          'ImagVibro' , ...         %6 
          'Null_1' , ...            %7
          'Null_2' , ...            %8
          'Att_1' , ...             %9
          'Att_2', ...              %10
          'Imag_vs_Att'};           %11
      
cvecs  = {[ 1  0  0  0  0  0  0  0  0  0  0 ], ...   % 1
          [ 0  1  0  0  0  0  0  0  0  0  0 ], ...   % 2
          [ 0  0  1  0  0  0  0  0  0  0  0 ], ...   % 3
          [ 0  0  0  1  0  0  0  0  0  0  0 ], ...   % 4
          [ 0  0  0  0  1  0  0  0  0  0  0 ], ...   % 5
          [ 0  0  0  0  0  1  0  0  0  0  0 ], ...   % 6 
          [ 0  0  0  0  0  0  1  0  0  0  0 ], ...   % 7
          [ 0  0  0  0  0  0  0  1  0  0  0 ], ...   % 8
          [ 0  0  0  0  0  0  0  0  1  0  0 ], ...   % 9
          [ 0  0  0  0  0  0  0  0  0  1  0 ], ...   %10
          [ 0  0  0  2  2  2  0  0 -3 -3  0 ]};      % 11
%         
      
del=1; % Delete existing contrasts (1=yes)
% were multiple regressors included in 1st level (step 1)?
n_hm=0;   % number of head motion parameters from realignment (step 4 in B0_preprocessing)
n_cc=0;   % number of CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
    % if >0, "zeros" will be appended in design matrix
    
%%
%# step 3: 2nd level glm, performs one-sample ttests over .con images

exclSJs = [3 9 11 15 25]; % subjects to be excluded from analysis 
outputfolder_2nd = 'otherMeanSearchlights'; % folder that will contain the created SPM file
dir_1st   = {'D4_r2agf4d_meanStim', 'D4_r2agf4d_meanImag'}; % Name of corresponding first Level Analysis

cnames_2nd = {'Stim', 'Imag'};
    
%           'StimPress' , ...         %1
%           'StimFlutt' , ...         %2
%           'StimVibro' , ...         %3
%           'ImagPress' , ...         %4
%           'ImagFlutt' , ...         %5
%           'ImagVibro' , ...         %6 
%           'Null_1' , ...            %7
%           'Null_2' , ...            %8
%           'Att' , ...               %9
%           'Response' , ...          %10
%           'StimPress_vs_N' , ...         %11
%           'StimFlutt_vs_N' , ...         %12
%           'StimVibro_vs_N' , ...         %13
%           'ImagPress_vs_N' , ...         %14
%           'ImagFlutt_vs_N' , ...         %15
%           'ImagVibro_vs_N' , ...         %16
%           'StimPress_vs_Att' , ...         %17 %%%%%%%%%%%%%%
%           'StimFlutt_vs_Att' , ...         %18
%           'StimVibro_vs_Att' , ...         %19
%           'ImagPress_vs_Att' , ...         %20
%           'ImagFlutt_vs_Att' , ...         %21
%           'ImagVibro_vs_Att' , ...         %22
%           'allStim_vs_N' , ...         %23
%           'allImag_vs_N' , ...         %24
%           'Att_vs_N' , ...         %25 %%%%%%%%%%%%%%%%%%%%
%           'allStim' , ...         %26
%           'allImag' , ...         %27
%           'allStim_vs_Att' , ...         %28
%           'allImag_vs_Att' , ...         %29
%           'Att_vs_allImag' , ...         %30
%           'allImag_vs_allStim', ...         %31
%           'allStim_vs_allImag'};             %32

%%
%# step 4: 2nd level glm, performs Flexible Factorial Design over .con images
con_images=1:10; % bei 2xN: 1:9,b 2xA: 1:10
exclSJs = [3 9 11 15 25]; % subjects to be excluded from analysis
outputfolder_2nd_b = '2nd_level_27_z_FlexFact'; % folder that will contain the created SPM file
dir_1st_b   = {'1st_level_27_z_ds8wragf4d'}; % Name of corresponding first Level Analysis

%#####################################################
%#################### INPUT end ######################
%#####################################################
%% step 1: 1st level glm, cycle over subjects
if ismember(1,analysis_switch)
    for sj = sbjs %1:length(SJs)
        display(SJs{sj})
        subj_dir = fullfile(src_dir, SJs{sj});
        try
            load([subj_dir filesep log_folder filesep logfile_prefix SJs{sj} '.mat'],'onsets');
            C1_glm_1stLevel(SJs, sj, subj_dir, outputfolder_1st, prefix_func, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc);
        catch 
            display('###################################################################')
            display(['################## ' SJs{sj} ', ERROR ###################'])
            display('###################################################################')
        end
    end
end
%% step 2: contrasts for 1st level, cycle over subjects
if ismember(2,analysis_switch)
    for sj = sbjs %1:length(SJs)
        display(SJs{sj})
        try
            display(SJs{sj})
            subj_dir = fullfile(src_dir, SJs{sj});
            C2_contrast_1stLevel(subj_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc);
        catch
            display('###################################################################')
            display(['################## ' SJs{sj} ', ERROR ###################'])
            display('###################################################################')
        end
    end
end
%% step 3: 2nd level glm, One Sample TTest
if ismember(3,analysis_switch)
    SJin=SJs;
    SJin(exclSJs)=[];
    C3_glm_2ndLevel_OneSampleTTest(src_dir, SJin, outputfolder_2nd, dir_1st, cnames_2nd);
end
%% step 4: 2nd level glm, Flexible Factorial Design
if ismember(4,analysis_switch)
    SJin=SJs;
    SJin(exclSJs)=[];
    C4_glm_2ndLevel_FlexFact(src_dir, SJin, outputfolder_2nd_b, dir_1st_b, con_images, runs);
end