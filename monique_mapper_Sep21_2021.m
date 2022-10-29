%WGGGG
%
%
% aim - to process data from metpro at the single subject level and
% populaiton level
% date - Jan 25, 2021
% author - saggar@stanford.edu
% updated on Sep 21, 2021 to clean the code

%% load data

cd metpro_processed/fmriprep-v1.5.9/fmriprep; 
clear; clc; close all;

% setting path for FSL commands to work
%setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
%setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');

stimDir = 'tda-metpro/data/StimTimes';
bx = readtable('tda-metpro/data/data_MetPro2_50subj_accuracy_anxiety_drugGuess.xlsx');
bx.STAI_median = (bx.STAI_delta_T3vsT1<5.1);

tasks = {'wm_run-2','wm_run-3'};

fd_thr = 0.2; 
skip_frames = 10; % initial TRs to skip

% for parcellation 
parcellation = 1;
parcellation_schemes = '';
data_parcel =  [];
conn = [];

%network based ciu
parcel_id = readtable('tda-metpro/data/schaefer_parcellation.xlsx');
network_names = {'Vis','SomMot','DAN','VAN','Lim','FP','DMN'};
runType = 'combined';
metricType = 'euclidean';

cmap = [
        47,101,165 % inst
        106,229,92 % rest
        251,227,86 % wm
        226,134,54 % vidoe
        174,52,40 %math
        ]./255;
cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
%% subject wise analysis - computing activation and Mapper statistics

networkActivation_inst_r1 = [];
networkActivation_safe_1bk_r1 = [];
networkActivation_safe_3bk_r1 = [];
networkActivation_thrt_1bk_r1 = [];
networkActivation_thrt_3bk_r1 = [];

networkActivation_inst_r2 = [];
networkActivation_safe_1bk_r2 = [];
networkActivation_safe_3bk_r2 = [];
networkActivation_thrt_1bk_r2 = [];
networkActivation_thrt_3bk_r2 = [];


for sub = 1:1:length(bx.id)
   sbj_name = sprintf('sub-%d', bx.id(sub));
   fprintf(1,'Processing %s \n', sbj_name);

   % get timing
   [timing_r(:,1), timing_r(:,2)] = readStimTiming(stimDir, bx.id(sub), skip_frames);

   % get data ready for mapper
   data_across_runs = [];
   timing_across_runs = [];
   runidx = [];
   numBadFrames_runs = [];
   for t = 1:1:length(tasks)
       fprintf(1,'Processing %s, task = %s\n', sbj_name, tasks{t});
       matVar = load(sprintf('%s/func/%s_task-%s_space-MNI152NLin6Asym_desc-preproc_bold_nuiReg_filtered_fd_02.mat', sbj_name, sbj_name, tasks{t}), 'matVar');
       matVar = matVar.matVar;
       
       parcel_cort = matVar.data_parcel_400;
       parcel_all = parcel_cort;
       data_across_runs = [data_across_runs; zscore(parcel_all(matVar.goodFrames,:))];
       timing_across_runs = [timing_across_runs; timing_r(matVar.goodFrames,t)]; 
       numBadFrames_runs = [numBadFrames_runs; length(matVar.badFrames)./(length(matVar.goodFrames)+length(matVar.badFrames))];
       runidx = [runidx; ones(length(matVar.goodFrames),1)*t];
   end
   
   networkData = [];
   for n = 1:1:length(unique(parcel_id.networkid))
       tmp = mean(data_across_runs(:,parcel_id.networkid==n),2); 
       networkData(:,n) = tmp;
       networkActivation_inst_r1(sub,n) = mean(networkData(timing_across_runs==1 & runidx == 1,n));
       networkActivation_safe_1bk_r1(sub,n) = mean(networkData(timing_across_runs==2 & runidx == 1,n));
       networkActivation_safe_3bk_r1(sub,n) = mean(networkData(timing_across_runs==3 & runidx == 1,n));
       networkActivation_thrt_1bk_r1(sub,n) = mean(networkData(timing_across_runs==4 & runidx == 1,n));
       networkActivation_thrt_3bk_r1(sub,n) = mean(networkData(timing_across_runs==5 & runidx == 1,n));

       networkActivation_inst_r2(sub,n) = mean(networkData(timing_across_runs==1 & runidx == 2,n));
       networkActivation_safe_1bk_r2(sub,n) = mean(networkData(timing_across_runs==2 & runidx == 2,n));
       networkActivation_safe_3bk_r2(sub,n) = mean(networkData(timing_across_runs==3 & runidx == 2,n));
       networkActivation_thrt_1bk_r2(sub,n) = mean(networkData(timing_across_runs==4 & runidx == 2,n));
       networkActivation_thrt_3bk_r2(sub,n) = mean(networkData(timing_across_runs==5 & runidx == 2,n));       
       
   end
  
   %run mapper
   try
    matVar_param_sweep = runBDLMapper_wrapper(data_across_runs, networkData, timing_across_runs, runidx, sbj_name, runType, metricType);    
    save(sprintf('tda-metpro/output/matVar_param_sweep_sub_%s_fd_02_%s', sbj_name, runType), 'matVar_param_sweep','data_across_runs','networkData','timing_across_runs','runidx');
   catch
    fprintf(2,'Error in main loop, ln 80 for sub=%s\n', sbj_name);
    continue;
   end

end

clc;
%% preparing data for mixed effects modeling
% Sep 20, 2022
networkActivation_safe_1bk = (networkActivation_safe_1bk_r1 + networkActivation_safe_1bk_r2)./2;
networkActivation_safe_3bk = (networkActivation_safe_3bk_r1 + networkActivation_safe_3bk_r2)./2;
networkActivation_thrt_1bk = (networkActivation_thrt_1bk_r1 + networkActivation_thrt_1bk_r2)./2;
networkActivation_thrt_3bk = (networkActivation_thrt_3bk_r1 + networkActivation_thrt_3bk_r2)./2;


%% model-driven analysis
mph_idx = 1:25;
pla_idx = 26:50;
dmn_net = 7;
frp_net = 6;
vis_net = 1;
som_net = 2;

%load based
alpha_r1_dmn = ((networkActivation_safe_3bk_r1(:,dmn_net)-networkActivation_safe_1bk_r1(:,dmn_net)) + (networkActivation_thrt_3bk_r1(:,dmn_net)-networkActivation_thrt_1bk_r1(:,dmn_net)))/2;
alpha_r2_dmn = ((networkActivation_safe_3bk_r2(:,dmn_net)-networkActivation_safe_1bk_r2(:,dmn_net)) + (networkActivation_thrt_3bk_r2(:,dmn_net)-networkActivation_thrt_1bk_r2(:,dmn_net)))/2;
beta_r1_dmn = (networkActivation_thrt_3bk_r1(:,dmn_net)+networkActivation_thrt_1bk_r1(:,dmn_net))/2 - (networkActivation_safe_3bk_r1(:,dmn_net)+networkActivation_safe_1bk_r1(:,dmn_net))/2;
beta_r2_dmn = (networkActivation_thrt_3bk_r2(:,dmn_net)+networkActivation_thrt_1bk_r2(:,dmn_net))/2 - (networkActivation_safe_3bk_r2(:,dmn_net)+networkActivation_safe_1bk_r2(:,dmn_net))/2;

alpha_r1_frp = ((networkActivation_safe_3bk_r1(:,frp_net)-networkActivation_safe_1bk_r1(:,frp_net)) + (networkActivation_thrt_3bk_r1(:,frp_net)-networkActivation_thrt_1bk_r1(:,frp_net)))/2;
alpha_r2_frp = ((networkActivation_safe_3bk_r2(:,frp_net)-networkActivation_safe_1bk_r2(:,frp_net)) + (networkActivation_thrt_3bk_r2(:,frp_net)-networkActivation_thrt_1bk_r2(:,frp_net)))/2;
beta_r1_frp = (networkActivation_thrt_3bk_r1(:,frp_net)+networkActivation_thrt_1bk_r1(:,frp_net))/2 - (networkActivation_safe_3bk_r1(:,frp_net)+networkActivation_safe_1bk_r1(:,frp_net))/2;
beta_r2_frp = (networkActivation_thrt_3bk_r2(:,frp_net)+networkActivation_thrt_1bk_r2(:,frp_net))/2 - (networkActivation_safe_3bk_r2(:,frp_net)+networkActivation_safe_1bk_r2(:,frp_net))/2;

alpha_r1_vis = ((networkActivation_safe_3bk_r1(:,vis_net)-networkActivation_safe_1bk_r1(:,vis_net)) + (networkActivation_thrt_3bk_r1(:,vis_net)-networkActivation_thrt_1bk_r1(:,vis_net)))/2;
alpha_r2_vis = ((networkActivation_safe_3bk_r2(:,vis_net)-networkActivation_safe_1bk_r2(:,vis_net)) + (networkActivation_thrt_3bk_r2(:,vis_net)-networkActivation_thrt_1bk_r2(:,vis_net)))/2;
beta_r1_vis = (networkActivation_thrt_3bk_r1(:,vis_net)+networkActivation_thrt_1bk_r1(:,vis_net))/2 - (networkActivation_safe_3bk_r1(:,vis_net)+networkActivation_safe_1bk_r1(:,vis_net))/2;
beta_r2_vis = (networkActivation_thrt_3bk_r2(:,vis_net)+networkActivation_thrt_1bk_r2(:,vis_net))/2 - (networkActivation_safe_3bk_r2(:,vis_net)+networkActivation_safe_1bk_r2(:,vis_net))/2;

alpha_r1_som = ((networkActivation_safe_3bk_r1(:,som_net)-networkActivation_safe_1bk_r1(:,som_net)) + (networkActivation_thrt_3bk_r1(:,som_net)-networkActivation_thrt_1bk_r1(:,som_net)))/2;
alpha_r2_som = ((networkActivation_safe_3bk_r2(:,som_net)-networkActivation_safe_1bk_r2(:,som_net)) + (networkActivation_thrt_3bk_r2(:,som_net)-networkActivation_thrt_1bk_r2(:,som_net)))/2;
beta_r1_som = (networkActivation_thrt_3bk_r1(:,som_net)+networkActivation_thrt_1bk_r1(:,som_net))/2 - (networkActivation_safe_3bk_r1(:,som_net)+networkActivation_safe_1bk_r1(:,som_net))/2;
beta_r2_som = (networkActivation_thrt_3bk_r2(:,som_net)+networkActivation_thrt_1bk_r2(:,som_net))/2 - (networkActivation_safe_3bk_r2(:,som_net)+networkActivation_safe_1bk_r2(:,som_net))/2;


alpha_dmn = (alpha_r1_dmn + alpha_r2_dmn)./2
beta_dmn = (beta_r1_dmn + beta_r2_dmn)./2

alpha_frp = (alpha_r1_frp + alpha_r2_frp)./2
beta_frp = (beta_r1_frp + beta_r2_frp)./2

alpha_vis = (alpha_r1_vis + alpha_r2_vis)./2
beta_vis = (beta_r1_vis + beta_r2_vis)./2

alpha_som = (alpha_r1_som + alpha_r2_som)./2
beta_som = (beta_r1_som + beta_r2_som)./2

%% all networks
%network_names = {'Vis','SomMot','DAN','VAN','Lim','FP','DMN'};
mph_idx = 1:25;
pla_idx = 26:50;

for net = 1:1:length(network_names)
    alpha_r1_net = ((networkActivation_safe_3bk_r1(:,net)-networkActivation_safe_1bk_r1(:,net)) + (networkActivation_thrt_3bk_r1(:,net)-networkActivation_thrt_1bk_r1(:,net)))/2;
    alpha_r2_net = ((networkActivation_safe_3bk_r2(:,net)-networkActivation_safe_1bk_r2(:,net)) + (networkActivation_thrt_3bk_r2(:,net)-networkActivation_thrt_1bk_r2(:,net)))/2;
    beta_r1_net = (networkActivation_thrt_3bk_r1(:,net)+networkActivation_thrt_1bk_r1(:,net))/2 - (networkActivation_safe_3bk_r1(:,net)+networkActivation_safe_1bk_r1(:,net))/2;
    beta_r2_net = (networkActivation_thrt_3bk_r2(:,net)+networkActivation_thrt_1bk_r2(:,net))/2 - (networkActivation_safe_3bk_r2(:,net)+networkActivation_safe_1bk_r2(:,net))/2;
    
    alpha_net(:,net) = (alpha_r1_net + alpha_r2_net)./2;
    beta_net(:,net) = (beta_r1_net + beta_r2_net)./2;
    
end
bx.alpha_net_vis = alpha_net(:,1);
bx.alpha_net_som = alpha_net(:,2);
bx.alpha_net_dan = alpha_net(:,3);
bx.alpha_net_van = alpha_net(:,4);
bx.alpha_net_lim = alpha_net(:,5);
bx.alpha_net_fp = alpha_net(:,6);
bx.alpha_net_dmn = alpha_net(:,7);

bx.beta_net_vis = beta_net(:,1);
bx.beta_net_som = beta_net(:,2);
bx.beta_net_dan = beta_net(:,3);
bx.beta_net_van = beta_net(:,4);
bx.beta_net_lim = beta_net(:,5);
bx.beta_net_fp = beta_net(:,6);
bx.beta_net_dmn = beta_net(:,7);

%% trying different bandwidth parameters
for bw = [.01,.05,.08,.2]
    figure; raincloud_plot(alpha_frp(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width', bw);
    raincloud_plot(alpha_frp(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width', bw);
    set(gcf,'color','w'); title(sprintf('alpha-frp-%f',bw)); axis([-.8 .8 -2 4]);
end

%%
figure; raincloud_plot(beta_dmn(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',[]);
raincloud_plot(beta_dmn(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',[]);
set(gcf,'color','w'); title('beta_dmn'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(alpha_dmn(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(alpha_dmn(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('alpha_dmn'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(beta_frp(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(beta_frp(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('beta_frp'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(alpha_frp(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(alpha_frp(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('alpha_frp'); axis([-.8 .8 -2 4]);


figure; raincloud_plot(beta_vis(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(beta_vis(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('beta_vis'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(alpha_vis(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(alpha_vis(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('alpha_vis'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(beta_som(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(beta_som(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('beta_som'); axis([-.8 .8 -2 4]);

figure; raincloud_plot(alpha_som(mph_idx), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
raincloud_plot(alpha_som(pla_idx), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
set(gcf,'color','w'); title('alpha_som'); axis([-.8 .8 -2 4]);




%% subject-wise analysis - analyzing mapper output
% redid analyss on Feb 28, 2022 

p_task = [];
z_task = [];
p_run = [];
z_run = [];
p_load = [];
z_load = [];
p_anx = [];
z_anx = [];
p_net = [];
z_net = [];
corr_nets = [];
p_load_noinst=[];
corr_across_params_r = [];
corr_across_params_p = [];
for param =  1:1:25
    param
    for s = 1:1:length(bx.id)
       sbj_name = sprintf('sub-%d', bx.id(s));
       if exist(sprintf('~/Dropbox/Articles_In_Preparation/tda-metpro/output/matVar_param_sweep_sub_%s_fd_02_%s.mat', sbj_name, runType),'file')     
        load(sprintf('~/Dropbox/Articles_In_Preparation/tda-metpro/output/matVar_param_sweep_sub_%s_fd_02_%s.mat', sbj_name, runType));
        try
            mapperout = matVar_param_sweep(param);
        catch
            fprintf(2,'error for %s \n', sbj_name);
            continue;
        end
                
        p_task(param,s) = mean(participation_coef(mapperout.nodeBynode, mapperout.nodeComm));
        z_task(param,s) = mean(module_degree_zscore(mapperout.nodeBynode, mapperout.nodeComm));

        p_run(param,s) = mean(participation_coef(mapperout.nodeBynode, mapperout.nodeComm_run));
        z_run(param,s) = mean(module_degree_zscore(mapperout.nodeBynode, mapperout.nodeComm_run));

         p_load(param,s) = mean(participation_coef(mapperout.nodeBynode, mapperout.nodeComm_load));
        tmp = participation_coef(mapperout.nodeBynode, mapperout.nodeComm_load);
        p_load_noinst(param,s) = mean(tmp(mapperout.nodeComm_load~=1));
        z_load(param,s) = mean(module_degree_zscore(mapperout.nodeBynode, mapperout.nodeComm_load));

        p_anx(param,s) = mean(participation_coef(mapperout.nodeBynode, mapperout.nodeComm_anx));
        z_anx(param,s) = mean(module_degree_zscore(mapperout.nodeBynode, mapperout.nodeComm_anx));
 
        p_net(param,s) = mean(participation_coef(mapperout.nodeBynode, mapperout.nodeComm_net));
        z_net(param,s) = mean(module_degree_zscore(mapperout.nodeBynode, mapperout.nodeComm_net));
         
       else
        continue;   
       end


    end

end
clc;


%% utility functions
function [timing_r1_f, timing_r2_f] = readStimTiming(stimDir, sub, trsDiscard)
    
    numTrs = 225;
    timing_r1_f = ones(numTrs,1); % as we know there must be 225 TRs for each run
    timing_r2_f = ones(numTrs,1);
    
    inst = load(sprintf('%s/%d/%d.instruction.stimtimes.1D', stimDir, sub, sub));
    safe1 = load(sprintf('%s/%d/%d.safeone.stimtimes.1D', stimDir, sub, sub));
    safe3 = load(sprintf('%s/%d/%d.safethree.stimtimes.1D', stimDir, sub, sub));    
    thrt1 = load(sprintf('%s/%d/%d.thrtone.stimtimes.1D', stimDir, sub, sub));
    thrt3 = load(sprintf('%s/%d/%d.thrtthree.stimtimes.1D', stimDir, sub, sub));
    shock = load(sprintf('%s/%d/%d.shock.stimtimes.1D', stimDir, sub, sub));
    
    timing = [inst'; safe1'; safe3'; thrt1'; thrt3';];
    indx = [ones(size(inst,2),1); ones(size(safe1,2),1)*2; ones(size(safe3,2),1)*3; ones(size(thrt1,2),1)*4; ones(size(thrt3,2),1)*5];
    
    timing_r1 = timing(:,1);
    indx_r1 = indx;
    timing_r2 = timing(:,2);
    indx_r2 = indx;
    
    [i,j]=sort(timing_r1);
    timing_r1 = round(timing_r1(j)/2)+1; %diving by Tr to get timing in units of Trs
    indx_r1 = indx_r1(j);
    
    [i,j]=sort(timing_r2);
    timing_r2 = round(timing_r2(j)/2)+1; %diving by Tr to get timing in units of Trs and addingn 1 for starting to count at 1
    indx_r2 = indx_r2(j);
    
    for tr = 1:1:length(timing_r1)-1
       timing_r1_f(timing_r1(tr):timing_r1(tr+1)-1) = indx_r1(tr);
    end
    timing_r1_f(timing_r1(end):timing_r1(end)+23) = indx_r1(end);
    
    for tr = 1:1:length(timing_r2)-1
       timing_r2_f(timing_r2(tr):timing_r2(tr+1)-1) = indx_r2(tr);
    end
    timing_r2_f(timing_r2(end):timing_r2(end)+23) = indx_r2(end);
    
    timing_r1_f(1:trsDiscard) = [];
    timing_r2_f(1:trsDiscard) = [];
end



function [matVar_param_sweep] = runBDLMapper_wrapper(parcelData, networkData, timing, run_idx, sbj_name, run, metricType)

    % run mapper
    data_z = parcelData;
    nclust = 10;
    num_bin_clusters = nclust;
    gain_val = 70;

    json_thrval = 0.5;
    json_saveOrNot = 1;

    % write json for load based piecharts    
    % timing: 1 ins, 2 1safe, 3, 3safe, 4 1threat, 5 3threat
    timing_load = timing;
    timing_load(timing_load==1) = 1;
    timing_load(timing_load==2) = 2;
    timing_load(timing_load==4) = 2;
    timing_load(timing_load==3) = 3;
    timing_load(timing_load==5) = 3;
    
   
    % write json for threat based piecharts
    % timing: 1 ins, 2 1safe, 3, 3safe, 4 1threat, 5 3threat
    timing_anx = timing;
    timing_anx(timing_anx==1) = 1;
    timing_anx(timing_anx==2) = 2;
    timing_anx(timing_anx==3) = 2;
    timing_anx(timing_anx==4) = 3;
    timing_anx(timing_anx==5) = 3;
    
    metaInfo.name = {'safe','threat'};
    metaInfo.timing = timing;
    metaInfo.run = run_idx;
    metaInfo.load = timing_load;
    metaInfo.anx = timing_anx;
    metaInfo.net = networkData;
    
    metaInfo_orig = metaInfo;
    itr = 1;
    k=16;
    for g = 66:2:74
        for r = 14:2:22
            try
                [nodeTpMat, nodeBynode, tpMat, filter, infTrs] = runBDLMapper(data_z, metricType, r, g, k, num_bin_clusters);
                json_output = sprintf('~/Dropbox/Articles_In_Preparation/tda-metpro/d3-nodepies/data/Schaefer_p400_fd02_%s_run_%s_metric_%s_res%d_gain%d_numk%d_3bin_nclust%d.json', sbj_name, run, metricType, r, g, k, nclust);
                metaInfo = metaInfo_orig;
                if ~isempty(infTrs)                   
                   metaInfo.timing(infTrs) = [];
                   metaInfo.run(infTrs) = [];
                   metaInfo.load(infTrs) = [];
                   metaInfo.anx(infTrs) = [];
                   metaInfo.net(infTrs,:) = [];
                end
                matVar = createJSON_NodePie_cme_enhanced(nodeTpMat, nodeBynode, metaInfo, json_output, json_thrval, json_saveOrNot);    
                matVar.metaInfo = metaInfo;
                matVar.nodeTpMat = nodeTpMat;
                matVar.nodeBynode = nodeBynode;
                matVar.tpMat = tpMat;
                matVar.filter = filter;
                matVar.res_val = r;
                matVar.gain_val = gain_val;
                matVar.num_k = k;
                matVar.num_bin_clusters = nclust;
                matVar.metricType = metricType;
                matVar.infTrs = infTrs;
                matVar_param_sweep(itr) = matVar;
                itr = itr + 1;
                
            catch
                matVar = struct();
                matVar_param_sweep(itr) = matVar;
                itr = itr + 1;
                continue;
            end
        end
        
    end

    
    
    
end
function [nodeTpMat, nodeBynode, tpMat, filter, infTrs] = runBDLMapper(data, metricType, res_val, gain_val, num_k, num_bin_clusters)

    X = data;


    resolution = [res_val res_val];
    gain = gain_val;
    
    fprintf(1,'Estimating distance matrix\n');
    tic
    distMat = estimateDistance(X, metricType);
    toc
    
    %--Manish Version--%
    fprintf(1,'Estimating knn graph\n');
    tic
    % create knn graph, estimate geodesic distances, embed using cmdscale and apply mapper
    [knnGraphTbl, knnGraph_dense_bin, knnGraph_dense_wtd, knnGraph_dense_bin_conn, knnGraph_dense_wtd_conn]= createPKNNG_bdl(distMat, num_k);

    knn_g_wtd = graph(knnGraph_dense_bin_conn);

    % estimate geodesic distances
    dist_geo_wtd = round(distances(knn_g_wtd,'Method','positive'));
    toc
    
    % check for inf values in the geodesic distance and take TRs out if any
    % Manish, Jan 28, 2021
    infTrs = find(isinf(dist_geo_wtd(:,1)));
    if ~isempty(infTrs)
       fprintf(2, 'Need to remove some TRs (%s) as their geodesic distances are inf\n', num2str(infTrs));
       X(infTrs, :) = [];

        fprintf(1,'Estimating distance matrix\n');
        tic
        distMat = estimateDistance(X, metricType);
        toc

        %--Manish Version--%
        fprintf(1,'Estimating knn graph\n');
        tic
        % create knn graph, estimate geodesic distances, embed using cmdscale and apply mapper
        [knnGraphTbl, knnGraph_dense_bin, knnGraph_dense_wtd, knnGraph_dense_bin_conn, knnGraph_dense_wtd_conn]= createPKNNG_bdl(distMat, num_k);

        knn_g_wtd = graph(knnGraph_dense_bin_conn);

        % estimate geodesic distances
        dist_geo_wtd = round(distances(knn_g_wtd,'Method','positive'));
        toc
    end
    
    fprintf(1,'Estimating embedding\n');
    tic
    
    % embed using cmdscale
    [y,e] = cmdscale(dist_geo_wtd);

    filter = [y(:,1), y(:,2)];
    toc
    
    fprintf(1,'Running mapper\n');
    tic
    
    
    %[adja, adja_pruned, pts_in_vertex, pts_in_vertex_pruned] = mapper2d_bdl_nonmetric(distMat, filter, resolution, gain, num_bin_clusters);
    [adja, adja_pruned, pts_in_vertex, pts_in_vertex_pruned] = mapper2d_bdl_hex_binning(distMat, filter, resolution, gain, num_bin_clusters, 3); %using triangulation
    %--End Manish Version--%
    toc
    
    fprintf(1,'Creating final output\n');
    tic
    % creating matrices for d3 visualization
    numNodes = length(pts_in_vertex_pruned);
    numTp = size(X,1);
    nodeTpMat = zeros(numNodes, numTp);
    for node = 1:1:numNodes
        tmp = pts_in_vertex_pruned{node};
        nodeTpMat(node, tmp) = 1;
    end

    nodeBynode = adja_pruned;
    tpMat = getMatTp_wtd(nodeBynode, nodeTpMat);
    fprintf(1,'Done\n');
    toc

end

function distMat = estimateDistance(X, metricType)
    distMat = squareform(pdist(X, metricType));
end

function matVar = createJSON_NodePie_cme_enhanced(nodeTpMat, nodeBynode, metaInfo, outpath, thrval, saveOrNot)

    runInfo = metaInfo.run;
    loadInfo = metaInfo.load;
    anxInfo = metaInfo.anx;
    netData = metaInfo.net;
    metaInfo = metaInfo.timing;
    
    
    nodeDeg = degrees_und(nodeBynode);
    
    numNodes = size(nodeTpMat,1);
    numTRs = size(nodeTpMat,2);
    thr = 1; %thrval; % 1 std. dev. above the mean for rsn metaInfor

    d = [];
    d.directed = 0;
    d.multigraph = 0;
    d.graph = {};
    tmp = sum(nodeTpMat,2);
    numLinks = 1;
    prop_nodes = [];
    prop_nodes_run = [];
    prop_nodes_load = [];
    prop_nodes_anx = [];
    prop_nodes_net = [];

    numMetaGroups = length(unique(metaInfo));
    numMetaGroups_run = length(unique(runInfo));
    numMetaGroups_load = length(unique(loadInfo));
    numMetaGroups_anx = length(unique(anxInfo));
    numMetaGroups_net = size(netData,2)+1; % last for NOTA
    
    nodeComm = [];
    nodeComm_run = [];
    nodeComm_load = [];
    nodeComm_anx = [];
    nodeComm_net = [];
    
    
    for node = 1:1:numNodes
        
        d.nodes{node}.row_count = tmp(node);
        trs = find(nodeTpMat(node,:));
        d.nodes{node}.trs = sprintf('[%s]',num2str(find(nodeTpMat(node,:))));

        prop = zeros(numMetaGroups,1); 
        prop_run = zeros(numMetaGroups_run,1); 
        prop_load = zeros(numMetaGroups_load,1); 
        prop_anx = zeros(numMetaGroups_anx,1); 
        prop_net = zeros(numMetaGroups_net,1); 
        
        for tr = trs
            prop(metaInfo(tr)) = prop(metaInfo(tr)) + 1;
            prop_run(runInfo(tr)) = prop_run(runInfo(tr)) + 1;            
            prop_load(loadInfo(tr)) = prop_load(loadInfo(tr)) + 1;            
            prop_anx(anxInfo(tr)) = prop_anx(anxInfo(tr)) + 1;   
%             if any(netData(tr,:)>thr)
%                 prop_net(netData(tr,:)>thr) = prop_net(netData(tr,:)>thr) + 1;
%             else
%                 prop_net(end) = prop_net(end) + 1;
%             end
            [~,j]=max(netData(tr,:));
            prop_net(j) = prop_net(j) + 1;
            
        end
                    
        [~,nodeComm(node)] = max(prop);
        prop_nodes(node,:) = prop;

        [~,nodeComm_run(node)] = max(prop_run);
        prop_nodes_run(node,:) = prop_run;

        [~,nodeComm_load(node)] = max(prop_load);
        prop_nodes_load(node,:) = prop_load;        
        
        [~,nodeComm_anx(node)] = max(prop_anx);
        prop_nodes_anx(node,:) = prop_anx;        

        [~,nodeComm_net(node)] = max(prop_net);
        prop_nodes_net(node,:) = prop_net;        
        
        for grp = 1:1:numMetaGroups_net
            d.nodes{node}.proportions_net{grp}.group = grp;
            d.nodes{node}.proportions_net{grp}.value = prop_net(grp);
            d.nodes{node}.proportions_net{grp}.row_count = tmp(node);  
            
            if grp < 6
                d.nodes{node}.proportions{grp}.group = grp;
                d.nodes{node}.proportions{grp}.value = prop(grp);
                d.nodes{node}.proportions{grp}.row_count = tmp(node);  
            else
                d.nodes{node}.proportions{grp}.group = grp;
                d.nodes{node}.proportions{grp}.value = 0;
                d.nodes{node}.proportions{grp}.row_count = tmp(node);  
                
            end

            if grp == 1
                d.nodes{node}.proportions_one{grp}.group = grp;
                d.nodes{node}.proportions_one{grp}.value = nodeDeg(node)+1;
                d.nodes{node}.proportions_one{grp}.row_count = tmp(node);
            else
                d.nodes{node}.proportions_one{grp}.group = grp;
                d.nodes{node}.proportions_one{grp}.value = 0;
                d.nodes{node}.proportions_one{grp}.row_count = tmp(node);                
                
            end
        end
        
        for grp = 1:1:numMetaGroups_net
            if grp < 3
                d.nodes{node}.proportions_run{grp}.group = grp;
                d.nodes{node}.proportions_run{grp}.value = prop_run(grp);
                d.nodes{node}.proportions_run{grp}.row_count = tmp(node);                             
            else
                d.nodes{node}.proportions_run{grp}.group = grp;
                d.nodes{node}.proportions_run{grp}.value = 0;
                d.nodes{node}.proportions_run{grp}.row_count = tmp(node);                                             
            end
        end
        
        for grp = 1:1:numMetaGroups_net
            
            if grp < 4
                d.nodes{node}.proportions_load{grp}.group = grp;
                d.nodes{node}.proportions_load{grp}.value = prop_load(grp);
                d.nodes{node}.proportions_load{grp}.row_count = tmp(node);     
            else
                d.nodes{node}.proportions_load{grp}.group = grp;
                d.nodes{node}.proportions_load{grp}.value = 0;
                d.nodes{node}.proportions_load{grp}.row_count = tmp(node);     
                
            end
        end
        
        for grp = 1:1:numMetaGroups_net
            if grp < 4
                d.nodes{node}.proportions_anx{grp}.group = grp;
                d.nodes{node}.proportions_anx{grp}.value = prop_anx(grp);
                d.nodes{node}.proportions_anx{grp}.row_count = tmp(node);     
            else
                d.nodes{node}.proportions_anx{grp}.group = grp;
                d.nodes{node}.proportions_anx{grp}.value = 0;
                d.nodes{node}.proportions_anx{grp}.row_count = tmp(node);     
                
            end
        end
        
       
        d.nodes{node}.id = node-1;

        links = find(nodeBynode(node,:));
        if ~isempty(links)        
            for l = 1:1:length(links)
                d.links{numLinks}.source = node - 1;
                d.links{numLinks}.target = links(l) - 1;
                numLinks = numLinks + 1;
            end

        end

    end
    
    if saveOrNot == 1 % only save when asked!
        savejson('',d,outpath);
    end
    
    matVar.nodeComm = nodeComm;
    matVar.nodeComm_run = nodeComm_run;
    matVar.nodeComm_load = nodeComm_load;
    matVar.nodeComm_anx = nodeComm_anx;
    matVar.nodeComm_net = nodeComm_net;
    
    matVar.prop_nodes = prop_nodes;
    matVar.prop_nodes_run = prop_nodes_run;
    matVar.prop_nodes_load = prop_nodes_load;
    matVar.prop_nodes_anx = prop_nodes_anx;
    matVar.prop_nodes_net = prop_nodes_net;
    

end


