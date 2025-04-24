%% Axon to Well Amplitude and Spike Relation Test for Theta

%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\delta\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
% re_fs=data.re_fs;
re_fs=1000;
data=data.filtered_data;
t_rec=300;
re_t=0:1/re_fs:t_rec-(1/re_fs);

%3.5 SD
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat")

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;
%% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

rng('default')
nPermutations=500;

shuffledLFPAngles=zeros(nPermutations,length(LFPAngles));
for nPerm=1:nPermutations
    randPermAngles=randperm(length(LFPAngles));
    shuffledLFPAngles(nPerm,:)=LFPAngles(randPermAngles);
end

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
end
logicalValidLFPs=zeros(1,length(re_t));
logicalValidLFPs(validLFPIndex)=1;

%% Compare phase of LFP large amplitudes to E10 spikes

%ok channels in FID 6
% E10, A8, E11, C11

spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\M9_spikes.mat");
spikes=spikes.index;
fs=25000;
t=0:1/fs:t_rec-(1/fs);


spikes=round(remap(spikes,1,length(t),1,length(re_t)));
logicalSpikes=zeros(1,length(re_t));
logicalSpikes(spikes)=1;

wellSpikeAngles=LFPAngles(logicalSpikes&logicalValidLFPs);
wellSpikeAngles=[wellSpikeAngles-360,wellSpikeAngles];

figure
histogram(wellSpikeAngles,-360:60:360)
xticks(-360:60:360)
xlabel("Axon Theta Angle")
ylabel("Soma Spike Count")
ax=gca;
ax.FontSize=16;


%% Test function
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==5 & well_spike_dyn.regi==3);
MI_tbl=cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,5,"G10",targetElecs,well_spike_dyn,40,thresh_mult);

%% scatter plot SPB vs Burst Length
scatter_BL_v_SPB(well_spike_dyn,MI_tbl(MI_tbl.pval~="NA",:),t,re_t,logicalValidLFPs)
%% Regression tests
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==5 & well_spike_dyn.regi==4);
sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,5,"G10",targetElecs,well_spike_dyn,20,thresh_mult,...
    "D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\")
