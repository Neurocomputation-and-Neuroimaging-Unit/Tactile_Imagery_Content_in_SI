% this skript does the actual decoding and stuff ordered in
% decoding_beta_shuffels.m

function[reference] = BetaShuffles(beta_dir, outputDir, rando, beta_selection)

% welcome to the subskript where the decodings are many and the references
% ain't pretty

% CREATE THE OUTPUT DIRECTORY IF IT DOES NOT EXIST YET
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% CONFIGURE THE DECODING TOOLBOX

clear cfg
cfg = decoding_defaults;
cfg.software = 'SPM12';

% Specify where the results should be saved
cfg.results.overwrite     = 1;
cfg.results.dir           = outputDir;

% DATA SCALING
cfg.scale.method          = 'z'; %z-score for all voxels over all samples (2 conditions x 6 runs =12)
cfg.scale.estimation      = 'all'; % only training, only test data, or all

% SEARCHLIGHT SPECIFICATIONS
cfg.analysis              = 'roi';   % or 'roi'
cfg.searchlight.unit      = 'voxels'; % or mm 
cfg.searchlight.radius    = 4; % 4 voxel or mm
cfg.searchlight.spherical = 0;  % only useful for mm
% The amount of information you want to have printed on the screen
% 0: no, 1: normal, output, 2: all)
cfg.verbose               = 1;  

% Method and model parameters 
cfg.decoding.method = 'classification';
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};

cfg.plot_selected_voxels  = 0;
cfg.plot_design           = 0;

% MASKs

mask_path = 'C:\Users\saraw\Desktop\Masks\Null_CONJ_CONJ\005';
cd(mask_path);
masks = dir('*.nii');
for ma=1:length(masks)
    cfg.files.mask(1,ma) = {masks(ma).name};
end

% DESIGN
labels = [-1 1];
cfg.files.label =  repmat(labels, 1, 6)';
cfg.files.chunk =  kron([1:6], ones(1,2))';
cfg.basic_checks.DoubleFilenameEntriesOk = 1;

% WE NEED:
% same decoding settings as with the results (IMPORTANT)
% one decoding per rando
% -> random selection of files for each run (times rando)
%    -> manual selection of betas instead of general (see other scripts)
%       -> predefined distribution of indexes of the betas
%           -> take a rando amount of these

if beta_selection == 0
    bins = [1:6; 11:16; 21:26; 31:36; 41:46; 51:56]; % zeile = run, spalte = regressor
elseif beta_selection == 1
% wenn nur stim
    bins = [1:3; 11:13; 21:23; 31:33; 41:43; 51:53]; % zeile = run, spalte = regressor
elseif beta_selection == 2
% wenn nur imag
    bins = [4:6; 14:16; 24:26; 34:36; 44:46; 54:56]; % zeile = run, spalte = regressor
elseif beta_selection == 3
    % imag and att -> 4 betas
    bins = [[4:6 9]; [14:16 19]; [24:26 29]; [34:36 39]; [44:46 49]; [54:56 59]];
end

reference = [];

for r = 1:rando
    % create random indices for betas
    ind_rand_betas = [];
    for n = 1:size(bins,1)
        ind_rand_betas = [ind_rand_betas randsample(bins(n,:),2)];
    end
        
    %Select Files by using a filter
    clear i
    f = {};
    for i = 1:length(ind_rand_betas)
        if ismember(i, [1 2])
            f{i} = [fullfile(beta_dir), ['\beta_000' num2str(ind_rand_betas(i)) '.nii']];
        elseif ismember(i, [3:12])
            f{i} = [fullfile(beta_dir), ['\beta_00' num2str(ind_rand_betas(i)) '.nii']];
        end
    end
    cfg.files.name = f';

    % define design
    cfg.design = make_design_cv(cfg);
    display_design(cfg);

    % Run decoding
    [references, cfg] = decoding(cfg);
    references = references(1,1).accuracy_minus_chance.output;

    % get references
    reference = [reference references];
end
