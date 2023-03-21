% this skript does the actual decoding and stuff ordered in
% decoding_beta_shuffels.m

function[reference] = BetaShuffles(beta_dir, outputDir, rando, beta_selection, labelnames)

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

cfg.results.output = {'accuracy','confusion_matrix'};

cfg.plot_selected_voxels  = 0;
cfg.plot_design           = 0;

mask_path = 'C:\Users\saraw\Desktop\newest_masks\GoodImag_model\controlClustersConjFWE005';
cd(mask_path);
masks = dir('*.nii');
for ma=1:length(masks)
    cfg.files.mask(1,ma) = {masks(ma).name};
end


% DESIGN
labels = [1 2 3 4 5 6];%[-1 1];
cfg.files.label =  repmat(labels, 1, 6)';
cfg.files.chunk =  kron([1:6], ones(1,length(labels)))';
cfg.basic_checks.DoubleFilenameEntriesOk = 1;

reference = [];



            load(fullfile(beta_dir, 'regressor_names.mat'));
            train1=find(strcmp(regressor_names(1,:),labelnames(1,1)));
            train2=find(strcmp(regressor_names(1,:),labelnames(1,2)));
            train3=find(strcmp(regressor_names(1,:),labelnames(1,3)));
                    train4=find(strcmp(regressor_names(1,:),labelnames(1,4)));
                    train5=find(strcmp(regressor_names(1,:),labelnames(1,5)));
                    train6=find(strcmp(regressor_names(1,:),labelnames(1,6)));

            f= {};
            filename='beta_0000';
            c       = []; % chunk
            l       = []; % label

            for i = 1:length(train1)
                temp=num2str(train1(i));
                f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                c(length(c)+1)  = i;
                l(length(l)+1)  = 1;


                temp=num2str(train2(i));
                f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                c(length(c)+1)  = i;
                l(length(l)+1)  = 2;


                temp=num2str(train3(i));
                f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                c(length(c)+1)  = i;
                l(length(l)+1)  = 3;


                            temp=num2str(train4(i));
                            f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                            c(length(c)+1)  = i;
                            l(length(l)+1)  = 4;

                
                            temp=num2str(train5(i));
                            f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                            c(length(c)+1)  = i;
                            l(length(l)+1)  = 5;

                
                            temp=num2str(train6(i));
                            f{length(f)+1} = spm_select('FPList',fullfile(beta_dir), [filename(1:end-length(temp)) temp '.(nii|img)']);
                            c(length(c)+1)  = i;
                            l(length(l)+1)  = 6;

            end

            cfg.files.name = f';
            cfg.files.chunk  = c';
            cfg.files.label  = l';


            % === Automatic Creation ===
            % This creates the leave-one-run-out cross validation design:

            % display_design(cfg);
            cfg.design = make_design_cv(cfg);


n_perms = rando;
combine = 0;   % see make_design_permutations how you can run all analysis in one go, might be faster but takes more memory
designs = make_design_permutation(cfg,n_perms,combine);

%% Run all permutations in a loop
% With small tricks to make it run faster (reusing design figure, loading 
% data once using passed_data), renaming the design figure, and to display
% the current permutation number in the title of the design figure)

cfg.fighandles.plot_design = figure(); % open one figure for all designs
passed_data = []; % avoid loading the same data multiple times by looping it

for i_perm = 1:n_perms

    dispv(1, 'Permutation %i/%i', i_perm, n_perms)
    
    cfg.design = designs{i_perm};
    cfg.results.filestart = ['perm' sprintf('%04d',i_perm)];
    
    set(cfg.fighandles.plot_design, 'name', sprintf('Permutation %i/%i', i_perm, n_perms)); % to know where we are
    if ~strcmp(cfg.analysis, 'searchlight') && i_perm > 1
        cfg.plot_selected_voxels = 0; % switch off after the first time, drawing takes some time
    end
    
    % do the decoding for this permutation
    [references, cfg, passed_data] = decoding(cfg, passed_data); % run permutation
    
    % rename  design figures to start with the current permutation number
    designfiles = dir(fullfile(cfg.results.dir, 'design.*'));
    for design_ind = 1:length(designfiles)
        movefile(fullfile(cfg.results.dir, designfiles(design_ind).name), ...
        fullfile(cfg.results.dir, [cfg.results.filestart '_' designfiles(design_ind).name]));
    end
% end

            %
            % %$ Run DECODING
            % results = decoding(cfg);
            %
            %     % define design
%             %     cfg.design = make_design_cv(cfg);
%             display_design(cfg);
% 
%             % Run decoding
%             [references, cfg] = decoding(cfg);
    references = references(1,1).accuracy.output;
    
    % get references
    reference = [reference references];
end
%     end
end
