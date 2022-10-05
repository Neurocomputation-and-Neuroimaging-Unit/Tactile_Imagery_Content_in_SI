% Evgeniya Kirilina
% FU Berlin, 29.01.2014
% Generating auxiliary regressors based on PCA analysis of CSF and WM
% signals

% A CompCor method as described in Behzadi et al., 2007
% A very brief description can also be found in Kirilina et al., 2015
% ('The quest for the best') and in Limanowski & Blankenburg (2015)
%--------------------------------------------------------------------------
% Kirilina2015:         8mm, 95% 95%
% Lemanowski2015:       4mm, 80% 90%

clear all; 
close all

data_dir     = 'D:\WICI\Data';

           
% SJs          = {'SJ01','SJ02_3','SJ03','SJ04','SJ05','SJ06','SJ07','SJ08',  'SJ09','SJ10',...
%                 'SJ11','SJ12',  'SJ13','SJ14','SJ15','SJ16','SJ17','SJ18_2','SJ19','SJ20',...
%                 'SJ21','SJ22',  'SJ23','SJ24'}; 

SJs = {'SJ36'};

%%
       
for sj = 1:length(SJs)
    
    % directory in which mask images and  the resulting *.mat files are written 
    DirR = fullfile(data_dir, SJs{sj},  'PCA_4mm8090'); 
    if ~exist(DirR, 'dir')
        mkdir(DirR)
    end
    cd(DirR)
    
    % directory path where subject-specific structural image is stored
    MaskPath    = fullfile(data_dir, SJs{sj},'T1');                          
    MaskCSF     = spm_select('FPList',MaskPath,'^c3sV.*\.(img|nii)$');  % select the segmented T1 CSF image
    MaskWM      = spm_select('FPList',MaskPath,'^c2sV.*\.(img|nii)$');  % select the segmented white matter image
     
     % first, smooth the CSF and WM structural segment.
     % Lemanowski used a 4 mm full width at half maximum Gaussian kernel.
     % Kirilina used a 8 mm FWHM Gaussian kernel

    fs = cellstr([MaskCSF; MaskWM]);     
    matlabbatch{1}.spm.spatial.smooth.data      = fs;
    matlabbatch{1}.spm.spatial.smooth.fwhm      = [4 4 4];
    matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
    matlabbatch{1}.spm.spatial.smooth.im        = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix    = 's4';
    
    spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);        
    clear matlabbatch fs

    cd(DirR)

% now select the smoothed CSF and WM mask for the next step
    SmoothedMaskCSF     = spm_select('FPList',MaskPath,'^s4c3sV.*\.(img|nii)$');  % select the segmented T1 CSF image
    SmoothedMaskWM      = spm_select('FPList',MaskPath,'^s4c2sV.*\.(img|nii)$');    


    for session=1:4
        % select the directory path in which functional EPIs are stored     
        EPIPath     = fullfile(data_dir, SJs{sj}, ['run0' num2str(session)]);    
        epi_filt    = '^raf.*\.(img|nii)$'; % Changed so use the slice time corrected, realigned EPIs!!!!!       
        Images      = cellstr(spm_select('FPList',EPIPath, epi_filt));     
                                                                           
        %create run-specific CSF and WM mask by multiplying the first image
        %of each run with the thresholding, smoothed CSF and WM segment   
        matlabbatch{1}.spm.util.imcalc.input = {
                                                Images{1,1}
                                                SmoothedMaskCSF
                                                };
        matlabbatch{1}.spm.util.imcalc.output = strcat(['run0' num2str(session) '_MaskCSF.img']);
        matlabbatch{1}.spm.util.imcalc.outdir = DirR;
        
        % set a  threschhold for CSF-Mask
        % Evgenia told me (YH) to use a really stringent threshold. in her
        % paper, she set the threshold to 0.95
        % YH had set the threshold to 0.8
        matlabbatch{1}.spm.util.imcalc.expression       = 'i1.*(i2>0.80)'; 
        matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask     = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
        
        spm_jobman('initcfg');
        spm_jobman('run', matlabbatch);
        clear matlabbatch        
        
        % the same game for white matter       
        matlabbatch{1}.spm.util.imcalc.input  = {
                                                Images{1,1}
                                                SmoothedMaskWM
                                                };
        matlabbatch{1}.spm.util.imcalc.output           = strcat(['run0' num2str(session) '_MaskWM.img']);
        matlabbatch{1}.spm.util.imcalc.outdir           = DirR;
        matlabbatch{1}.spm.util.imcalc.expression       = 'i1.*(i2>0.90)'; 
        matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask     = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
        
        spm_jobman('initcfg');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
        
        % get the header information for the images
        EPIs =  spm_vol(spm_select('FPList', EPIPath, epi_filt));
        mCSF =  spm_vol(spm_select('FPList', DirR, ['^run0' num2str(session) '_MaskCSF.*\.(img|nii)$']));
        mWM =   spm_vol(spm_select('FPList', DirR, ['^run0' num2str(session) '_MaskWM.*\.(img|nii)$']));
        
        % read in the image data
        Y(:,:,:,:)   =  spm_read_vols(EPIs);            % the value of each voxel across all acquisition points
        mCSFa(:,:,:) =  spm_read_vols(mCSF);            % same for the CSF mask
        mWMa(:,:,:)  =  spm_read_vols(mWM);             % same for the WM mask
        
        ss           = size(Y);                         % [x y z N volumes]
        
        % The next 3 lines had been commented out by YH or EK
        % Either 0 out empty voxels for CSF and WM together or separately
        % I (LV) don't what what difference is between the two
%         yCSFWM                                      =  reshape(Y, ss(1)*ss(2)*ss(3), ss(4));
%         yCSFWM((mCSFa == 0)&(mWMa == 0),:)          =  [];   % yCSFWM(find((mCSFa == 0)&(mWMa == 0)),:) =  [];   % discard voxels which neither belongs to CSF nor to WM
% 
%         % substract from each voxel its the time series mean  
%         yCSFWM                                      = yCSFWM - repmat( mean(yCSFWM,2), [1, size(yCSFWM,2)]); 
        
        yCSF                = reshape(Y, ss(1)*ss(2)*ss(3), ss(4));
        yCSF(mCSFa == 0,:)  = [];                       %  yCSF(find(mCSFa == 0),:) = [];
        yCSF                = yCSF-repmat( mean(yCSF,2), [1, size(yCSF,2)]);        
        
        yWM                 = reshape(Y, ss(1)*ss(2)*ss(3), ss(4));
        yWM(mWMa == 0,:)    = [];                       % yWM(find(mWMa == 0),:) = [];
        yWM                 = yWM-repmat( mean(yWM,2), [1, size(yWM,2)]);
         
        clear EPIPath EPIs Images 
        
        % Extraxt ROI-Data from CSF and WM
        
        % yCSF        = spm_get_data(Images,xSPM.XYZ(:,Q));
        % yWM         = spm_get_data(Images,xSPM.XYZ(:,Q));        
        % y           = spm_filter(SPM.xX.K,SPM.xX.W*y); filtering
        
        % Compute regional CSF response in terms of first five PCA
        %--------------------------------------------------------------------------        
        coeff           = pca(yCSF'*yCSF);
        signalCSF       = coeff(:,1:5);
        
        clear coeff
        
        R               = signalCSF;
        filenameCSF     = fullfile(DirR, ['run0' num2str(session) '_csf.mat']);
        save (filenameCSF, 'R');
        
        clear R
        
        %- Compute regional WM response in terms of first five PCA
        %--------------------------------------------------------------------------    
        coeff          = pca(yWM'*yWM);
        signalWM       = coeff(:,1:5);
        
        clear coeff
        
        R              = signalWM;
        filenameWM     = fullfile(DirR, ['run0' num2str(session) '_wm.mat']);
        save (filenameWM, 'R');
        
        clear R
         
        %cd(data_dir);
        
        %R =[signalCSF signalWM];     
                
        clear Y yCSF yWM;
        clear mCSFa mCSF signalCSF;
        clear mWMa mWM signalWM;
    end;
end;

%%
%-Set in structure
%--------------------------------------------------------------------------
% xY.y    = y;
% xY.u    = YCSF;
% xY.v    = v;
% xY.s    = s;
% 
% save(fullfile(SPM.swd,str),'YCSF','YWM')