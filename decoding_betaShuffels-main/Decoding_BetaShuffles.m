% this skript starts the decoding of data selected at random
% there will be n decodings with a differnt set of random data each
% the data used will be the 6 "main" regressors, one picked at random for
% each run
%
% afterwards REAL results will be loaded
% there will be a p-value assigned to the results which describes the
% chance that such a result occurs when using random data
%
% watch the runtime -> open more than one script to run some 
% simultaneously -> watch the CPU

%% BASICS

switchy = [3]; % 1: rando-decoding, 2: loading, 3: stats

src_dir      = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data';

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

sbjs = [1 2 4 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32]; 
% sbjs = [1 2 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32]; 

betaDir = '1st_level_D0_wragf4d';

outputDir = 'DX_wragf4d_SHUFFLE_Null_CONJ_CONJ005';

mask_path = 'C:\Users\saraw\Desktop\Masks\Null_CONJ_CONJ\005';
cd(mask_path);
mask_names = dir('*.nii');
masks = {};
for ma=1:length(mask_names)
    masks{ma} = [mask_path mask_names(ma).name];
end
roi_names = masks;

labelpairs = {'StimPress','StimFlutt'; ...
    'StimPress','StimVibro'; ...
    'StimFlutt','StimVibro'; ...
    'ImagPress','ImagFlutt'; ...
    'ImagPress','ImagVibro'; ...
    'ImagFlutt','ImagVibro'}; 

rando = 1000; % number of planned decodings on betas picked at random
seed = '1000_6_betas_Null_CONJ_CONJ005'; % so you dont overwrite your first references(if you have more than one skript running)
beta_selection = 0;% do we want the 6 (0), stim (1), imag (2) or Imag_Att (3)??

%% DECODING
% WE WANT: 
% one matrix of references (sj, roi, rando) in the src_dir
% one matrix of references per sj (roi, rando) in the outputDir
% WE NEED:
% a loop over SJs
% a decoding per ROI (we can handle that in BetaShuffles)
% a loop over the randos (handle in subscript? -> less defining neccessary)
% a thing that stores all references

if ismember(1, switchy)
    reference = zeros(length(sbjs), length(roi_names), rando);
    
    for sj = sbjs
    
        s = find(sbjs == sj);
    
        display(['Step 4, Decoding: ' SJs{sj} '; ' num2str(rando) ' random iterations'])
        beta_path    = fullfile(src_dir, SJs{sj}, betaDir);
        output_path  = fullfile(src_dir, SJs{sj}, outputDir);
    
        this_reference = BetaShuffles(beta_path, output_path, rando, beta_selection); 
    
        writematrix(this_reference, [output_path '\references_ROIs_x_rando_' seed]);
    
        reference(s,:,:) = this_reference;
    
    end
    
    save([src_dir '\references_ROIs_x_rando_' seed '.mat'], 'reference');

end

%% LOADING
% take results per person & roi
% take references per person & roi (& rando)

outputDir4 = 'D4_wragf4d_Null_CONJ_CONJ005';

if ismember(2, switchy)

    results_temp = [];
    results = zeros(length(sbjs),length(roi_names),6); %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Attention: actualize to your seeds!!!
    reference_file = {'references_ROIs_x_rando_1000_6_betas_Null_CONJ_CONJ005', 'references_ROIs_x_rando_a1000permsStim', 'references_ROIs_x_rando_a1000permsImag', 'references_ROIs_x_rando_1000_Imag_Att'};
    task = 0; % 0: just take the first try, 1: the only-stim-file, 2: the only-imag-file, 3: the 5000er file and all are tested against it 4: xClass on 5000er
    imag_att = [];
    xClass_Stim_Imag = [];
    xClass_Imag_Stim = [];

    % get results and references
    for sj = sbjs
        s = find(sbjs == sj); % you know the struggle

        for lap = 1:size(labelpairs, 1)
        
            %results
%             temp_results = load([src_dir '\' SJs{sj} '\' outputDir4 '_' labelpairs{lap, 1} '_' labelpairs{lap, 2} '\plain_results.mat']);
%             results_temp = [results_temp; temp_results(1,1).saveable'];
%             results(s,:,lap) = temp_results(1,1).saveable'; % (sj, roi, labelpair)

        end

        %att
%         folder_names_att = {'\D4_wragf4d_ImagPress_Att', '\D4_wragf4d_ImagFlutt_Att', '\D4_wragf4d_ImagVibro_Att'};
%         outputDir_att = {[src_dir '\' SJs{sj} '\D4_wragf4d_Imag_Att']};
%         temp_imag_att = averageResults(src_dir, sj, folder_names_att, outputDir_att);
%         imag_att = [imag_att; temp_imag_att'];

        %xClass
%         folder_names_xClass_Stim_Imag = {'\D5_XClass_wragf4d_TRAIN_StimPress_StimFlutt_TEST_ImagPress_ImagFlutt', ...
%                                          '\D5_XClass_wragf4d_TRAIN_StimPress_StimVibro_TEST_ImagPress_ImagVibro', ...
%                                          '\D5_XClass_wragf4d_TRAIN_StimFlutt_StimVibro_TEST_ImagFlutt_ImagVibro'};
%         folder_names_xClass_Imag_Stim = {'\D5_XClass_wragf4d_TRAIN_ImagPress_ImagFlutt_TEST_StimPress_StimFlutt', ...
%                                          '\D5_XClass_wragf4d_TRAIN_ImagPress_ImagVibro_TEST_StimPress_StimVibro', ...
%                                          '\D5_XClass_wragf4d_TRAIN_ImagFlutt_ImagVibro_TEST_StimFlutt_StimVibro'};
%         outputDir_x = {[src_dir '\' SJs{sj} '\D5_wragf4d_xClass_Stim_Imag'], [src_dir '\' SJs{sj} '\D5_wragf4d_xClass_Imag_Stim']};
%         xClass_Stim_Imag = [xClass_Stim_Imag; averageResults(src_dir, sj, folder_names_xClass_Stim_Imag, outputDir_x(1))'];
%         xClass_Imag_Stim = [xClass_Imag_Stim; averageResults(src_dir, sj, folder_names_xClass_Imag_Stim, outputDir_x(2))'];

%             results = data3D;


        %references
        if task == 0
            reference = load([src_dir '\' cell2mat(reference_file(1))]);
            reference = reference(1,1).reference;
        elseif task == 1
            reference = load([src_dir '\' cell2mat(reference_file(2))]);
            reference = reference(1,1).reference;
        elseif task == 2
            reference = load([src_dir '\' cell2mat(reference_file(3))]);
            reference = reference(1,1).reference;
        elseif ismember(task, [3 4])
            reference = load([src_dir '\' cell2mat(reference_file(4))]);
            reference = reference(1,1).reference;
        end
    end
end

%% STATS  
% give p-value which describes probability that result could have been
% created by the random procedure 
% e.g.: 100 randos -> minimal p = 0.01, 1000 randos -> minimal p = 0.001
task = 0;

if ismember(3, switchy)

    results = data3D;
%     results(12, :, :) = [];
%     reference(12, :, :) = [];
%     task = 3
%     reference(3, :, :) = [];

    rando = rando;

    p = zeros(length(roi_names), 1);
    pA = p;
    pB = p;
    pC = p;
    p1 = p;
    p2 = p;
    p3 = p;
    p4 = p;

    % do that stats-stuff
    for roi = 1:size(results,2)

        current_results = reshape(results(:,roi,:),[length(sbjs) 6]);
        mean_results_stim = mean(current_results(:,[1:3]),'all');
        mean_results_imag = mean(current_results(:,[4:6]),'all');
        imag1 = mean(current_results(:,[4]),'all');
        imag2 = mean(current_results(:,[5]),'all');
        imag3 = mean(current_results(:,[6]),'all');
        
     

%                 mean_imag_att = mean(results(:,roi), 'all');
%         mean_xClass_stimImag = mean(xClass_Stim_Imag(:,roi), 'all');
%         mean_xClass_imagStim = mean(xClass_Imag_Stim(:,roi), 'all');

        current_reference = mean(reshape(reference(:,roi,:),[length(sbjs) rando]),1);
        sz_ref = size(current_reference);

% Result should be p(roi) with means for each category stim or imag /
% imag-att
        if task == 0
            % the main plot
            p1(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_results_stim,current_reference),2)+1);
            p2(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_results_imag,current_reference),2)+1);
%             p3(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_xClass_stimImag,current_reference),2)+1);
%             p4(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_xClass_imagStim,current_reference),2)+1);
            pA(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,imag1,current_reference),2)+1);
            pB(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,imag2,current_reference),2)+1);
            pC(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,imag3,current_reference),2)+1);

        elseif task == 1
            p(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_results_stim,current_reference),2)+1);
        elseif task == 2
            p(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_results_imag,current_reference),2)+1);
        elseif task == 3
            p3(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_imag_att,current_reference),2)+1);
        elseif task == 4
            p1(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_xClass_stimImag,current_reference),2)+1);
            p2(roi) = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,mean_xClass_imagStim,current_reference),2)+1);
        end
    end
end


