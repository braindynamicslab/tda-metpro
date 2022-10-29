%WGGGG
%
%
% aim - postprocess rsfMRI data after processing through fMRIprep
% author - saggar@stanford.edu
% date - 8.24.2020
% version - v1 (Aug 2020)
% Sep 14, 2020: updated to include parcellation as optional
% Oct 7, 2020: updates to include multiple parcellations and sub-cortical
% parcellation


%% data
cd metpro_processed/fmriprep-v1.5.9/fmriprep; 
clear; clc; close all;

% setting path for FSL commands to work
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');


tmp = dir('sub-*.html');
subjs = strrep({tmp(:).name},'.html',''); %subject names

tasks = {'wm_run-2','wm_run-3'};

fd_thr = 0.2; % for clinical data
skip_frames = 10; % initial TRs to skip

% for parcellation 
parcellation = 1;
parcellation_schemes = '';
data_parcel =  [];
conn = [];


% temporal filtering setup
TR = 2;
Fs = 1/TR;
Ny = Fs/2;
BW_Hz = [0.009 0.08];
bp_order = 4; % for a second order filter (2N for bandpass)
BW_N = BW_Hz/Ny; % dividing by Ny frequency as butter requires normalized Wn
[b_filter, a_filter] = butter(bp_order, BW_N);

%% subject wise analysis
for sub = 1:1:length(subjs)
   sbj_name = subjs{sub};
   fprintf(1,'Processing %s \n', sbj_name);
   
   for t = 1:1:length(tasks)
       fprintf(1,'Processing %s, task = %s\n', sbj_name, tasks{t});
       
       % confound file
       tmp = sprintf('%s/func/%s_task-%s_desc-confounds_regressors', sbj_name, sbj_name, tasks{t});
       cmd = sprintf('cp %s.tsv %s.txt', tmp, tmp);
       system(cmd);

       confounds = readtable([tmp '.txt'], 'TreatAsEmpty', 'n/a');
       
       % skipping initial frames
       confounds(1:skip_frames,:) = [];
       
       % scrubbing, also tagging one back and two forward frames
       fd = confounds.framewise_displacement;
       badFrames = find(fd > fd_thr);
       badFrames = unique([badFrames; badFrames-1; badFrames+1; badFrames+2]);
       
       badFramesZero = find(badFrames==0);
       badFrames(badFramesZero) = [];
       badFramesHigh = find(badFrames > size(confounds,1));
       badFrames(badFramesHigh) = [];       
             
       goodFrames = setdiff(1:size(confounds,1), badFrames);
       
       % got 6P regressors from fmriprep
       nuisance_reg = [confounds.trans_x, confounds.trans_y, confounds.trans_z, confounds.rot_x, confounds.rot_y, confounds.rot_z];
       
       % replacing NaN on first row with zeros
       nuisance_reg(isnan(nuisance_reg)) = 0; % making sure first value of derivatives is zero
       
       % put high motion frames as nan
       nuisance_reg(badFrames,:) = nan;
       
       % demean and detrend nuisance regressors
       nuisance_reg_colmean = nanmean(nuisance_reg); %for demeaning only take low motion frames
       nuisance_reg_demean = nuisance_reg - nuisance_reg_colmean; %but demean all frames
       
       % detrending using only low motion frames
       nuisance_reg_demean_detrend = detrend(nuisance_reg_demean, 'omitnan'); %detrend only low motion frames
       
       % read data 
       data = niftiread(sprintf('%s/func/%s_task-%s_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz', sbj_name, sbj_name, tasks{t}));
       
       % skipping initial data frames
       data(:,:,:,1:skip_frames) = [];
       
       data_info = niftiinfo(sprintf('%s/func/%s_task-%s_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz', sbj_name, sbj_name, tasks{t}));
       
       % creating a mask
       mask_filename = sprintf('%s/anat/%s_space-MNI152NLin6Asym_desc-brain_mask_2mm_bin.nii.gz', sbj_name, sbj_name);
       if ~exist(mask_filename, 'file')
        cmd = sprintf('fslmaths %s/anat/%s_space-MNI152NLin6Asym_desc-brain_mask.nii.gz -subsamp2 -bin %s/anat/%s_space-MNI152NLin6Asym_desc-brain_mask_2mm_bin', sbj_name, sbj_name, sbj_name, sbj_name);
        system(cmd);
       end       
       data_mask = niftiread(mask_filename);
              
       numX = size(data,1);
       numY = size(data,2);
       numZ = size(data,3);
       numT = size(data,4);
       
       data = reshape(data, [numX*numY*numZ, numT]);
       data = data'; %numT x voxels 
       
       data_mask = reshape(data_mask,  [numX*numY*numZ,1]);
       data_nz = double(data(:, data_mask==1)); %Only taking brain masked voxels       
       
       % putting censored frames out so they don't affect demean and
       % detrend       
       data_nz(badFrames,:) = nan;
       
       % demean data, column-wise demean
       data_nz = data_nz - nanmean(data_nz);
       
       % detrend, column-wise, for low motion frames, takes about 3 min 
       tic
       data_nz = detrend(data_nz, 'omitnan');
       toc
       
       % regress nuisance from the BOLD data
       R = zeros(size(data_nz));
       tic
       for vox = 1:1:size(data_nz,2)
        % only storing residuals
        [~,~,R(:,vox)] = regress(data_nz(:,vox), nuisance_reg_demean_detrend);        
       end
       toc
       
       % interpolation 
       R(badFrames,:) = interp1(goodFrames, R(goodFrames,:), badFrames);       
       mean_goodFrames = mean(R(goodFrames,:));
       
       temp = find(isnan(mean(R(badFrames,:),2)));
       R(badFrames(temp),:) = repmat(mean_goodFrames, numel(temp), 1);
       
       % temporal filtering
       R_filtered = filtfilt(b_filter, a_filter, R);
       
       if parcellation == 1
           % parcellation 1
           lbl = niftiread('Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
           lbl(data_mask==0)=[];
           data_parcel_100 = [];
           for pr = 1:1:max(lbl)
               data_parcel_100(:,pr) = mean(R_filtered(:, lbl==pr), 2);
           end
           
           lbl = niftiread('Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
           lbl(data_mask==0)=[];
           data_parcel_200 = [];
           for pr = 1:1:max(lbl)
               data_parcel_200(:,pr) = mean(R_filtered(:, lbl==pr), 2);
           end

           lbl = niftiread('Schaefer2018_300Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
           lbl(data_mask==0)=[];
           data_parcel_300 = [];
           for pr = 1:1:max(lbl)
               data_parcel_300(:,pr) = mean(R_filtered(:, lbl==pr), 2);
           end

           lbl = niftiread('Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
           lbl(data_mask==0)=[];
           data_parcel_400 = [];
           for pr = 1:1:max(lbl)
               data_parcel_400(:,pr) = mean(R_filtered(:, lbl==pr), 2);
           end
           
           lbl = niftiread('HarvardOxford-sub-maxprob-thr50-2mm.nii.gz');
           lbl(data_mask==0)=[];
           data_parcel_sc = [];
           for pr = 1:1:max(lbl)
               data_parcel_sc(:,pr) = mean(R_filtered(:, lbl==pr), 2);
           end
           
           % connectome
           tmp = data_parcel_400;
           tmp(badFrames,:) = [];
           conn = corr(tmp);           
       end

       
       
       % write nifti       
       data_to_write = zeros(numT, [numX*numY*numZ]);
       data_to_write(:, data_mask==1) = R_filtered;       
       data_to_write = data_to_write'; % voxels x numT
       data_to_write = reshape(data_to_write, [numX, numY, numZ, numT]);
       
       % changing data info
       data_info.Datatype = 'double';
       data_info.Filemoddate = datetime('now');
       data_info.ImageSize = size(data_to_write);
       data_info.Description = sprintf('data:nuisance_regressed:filtered:withbadFrames');
       niftiwrite(data_to_write, sprintf('%s/func/%s_task-%s_space-MNI152NLin6Asym_desc-preproc_bold_nuiReg_filtered_fd_02.nii.gz', sbj_name, sbj_name, tasks{t}), data_info, 'Compressed', 1);

       % write parcels and metainformation
       matVar.fd = fd;
       matVar.badFrames = badFrames;
       matVar.goodFrames = goodFrames;
       matVar.nuisance_reg = nuisance_reg;
       matVar.nuisance_reg_demean_detrend = nuisance_reg_demean_detrend;
       matVar.data_parcel_100 = data_parcel_100;
       matVar.data_parcel_200 = data_parcel_200;
       matVar.data_parcel_300 = data_parcel_300;
       matVar.data_parcel_400 = data_parcel_400;
       matVar.data_parcel_sc = data_parcel_sc;
       matVar.conn = conn;
       matVar.fd_thr = fd_thr;
       matVar.skip_frames = skip_frames;
       matVar.bp_order = bp_order;
       matVar.BW_Hz = BW_Hz;       
       save(sprintf('%s/func/%s_task-%s_space-MNI152NLin6Asym_desc-preproc_bold_nuiReg_filtered_fd_02.mat', sbj_name, sbj_name, tasks{t}), 'matVar');
        
       % write connectivity file
       save(sprintf('%s/func/%s_task-%s_funcConnectivity_fd_02.txt', sbj_name, sbj_name, tasks{t}), 'conn', '-ascii');

   end

    
end
