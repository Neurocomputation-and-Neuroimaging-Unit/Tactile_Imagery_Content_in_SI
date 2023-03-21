function[parameterC] = D4_Decoding_batch(beta_dir, output_dir, labelnames)
% This Batch Script first specifies what features should be decoded and
% executes the Decoding analysis.
% The resulting accuracy maps should later be normalized and smoothed

% CREATE THE OUTPUT DIRECTORY IF IT DOES NOT EXIST YET
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% CONFIGURE THE DECODING TOOLBOX

clear cfg
cfg = decoding_defaults;
cfg.software = 'SPM12';

% Specify where the results should be saved
cfg.results.overwrite     = 1;
cfg.results.dir           = output_dir;

% DATA SCALING
cfg.scale.method          = 'z'; %z-score for all voxels over all samples (2 conditions x 6 runs =12)
cfg.scale.estimation      = 'all'; % only training, only test data, or all

% SEARCHLIGHT SPECIFICATIONS
cfg.analysis              = 'roi';   % or 'roi'
cfg.searchlight.unit      = 'voxels'; % or mm 
cfg.searchlight.radius    = 5; % 4 voxel or mm
cfg.searchlight.spherical = 0;  % only useful for mm
% The amount of information you want to have printed on the screen
% 0: no, 1: normal, output, 2: all)
cfg.verbose               = 1;  

% Method and model parameters 
cfg.decoding.method = 'classification';
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; %-s 0 -t 0 -c 1 -b 0 -q

% OUTPUTS SPECIFICATION
cfg.results.output = {'accuracy','confusion_matrix'};%{'accuracy_minus_chance'};

% DISPLAY:
cfg.plot_selected_voxels  = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% This is by default set to 1, but if you repeat the same design again and again, it can get annoying...
cfg.plot_design           = 0;

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")

% cfg.files.mask = fullfile(beta_dir, 'mask.nii');

mask_path = 'C:\Users\saraw\Desktop\newest_masks\GoodImag_model\controlClustersConjUncorr_k200';
cd(mask_path);
masks = dir('*.nii');
for ma=1:length(masks)
    cfg.files.mask(1,ma) = {masks(ma).name};
end

%%% parameter selection

% cfg.parameter_selection.method = 'grid';
%   % grid search (currently the only implemented method)
% cfg.parameter_selection.parameters = {'-c'};
% cfg.parameter_selection.parameter_range = {[0.0001 0.001 0.01 0.1 1 10 100 1000]};

% % The following parameters are set as defaults, i.e., they don’t need to be called explicitly
% cfg.parameter_selection.format.name =  'string_number';
% cfg.parameter_selection.format.separator = ' ';
% cfg.parameter_selection.optimization_criterion = 'max'; % pick value with maximal accuracy

%%% feature transformation & selection
% 
% cfg.feature_transformation.method = 'PCA';
% cfg.feature_transformation.estimation = 'all';
% cfg.feature_transformation.critical_value = 0.3; % only keep components that explain at least 10 percent variance

% cfg.feature_selection.method = 'filter';
% cfg.feature_selection.filter = 'f0';
% cfg.feature_selection.n_vox = 'automatic';
   % nested CV to determine optimal number of features

%%% design

% Decoding DESIGN
labels = [1; 2; 3]; %-1; 1 %%%%%% ; 4; 5; 6

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% Extract all information for the cfg.files structure
cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir);

% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg);
display_design(cfg);

% Run decoding
[results] = decoding(cfg);

% cfg.decoding.train.classification.model_parameters
