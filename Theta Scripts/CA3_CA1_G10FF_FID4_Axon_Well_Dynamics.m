%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\Theta\4x 24574 210715 21div 210806_1_mat_files\G10.mat");
re_fs=data.re_fs;
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
nsamples_combine_thresh=(1/10)*re_fs; %1 cycle of fastest theta
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

%% Compare phase of axon LFP large amplutudes to E10 burst start
well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==5 & well_spike_dyn.channel_name=="M9"};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;

figure
wellBurstStartAngles=LFPAngles(logicalBurstStarts&logicalValidLFPs);
phaseBurstStartProb=histogram(wellBurstStartAngles,0:20:360,"Normalization","probability");
MI_BurstStart=spike_amplitude_MI(phaseBurstStartProb.Values);
wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];

figure
histogram(wellBurstStartAngles,-360:20:360)
xticks(-360:60:360)
xlabel("Axon Theta Angle")
ylabel("Soma Burst Start Count")
ax=gca;
ax.FontSize=16;

figure
histogram(wellBurstStartAngles,-360:20:360,"Normalization","probability")
xticks(-360:60:360)
xlabel("Axon Theta Angle")
ylabel("Soma Burst Start Probability")
ax=gca;
ax.FontSize=16;

% plots amplitude in axon against target burst start
wellBurstStartAmp=LFPAmplitude(logicalBurstStarts&logicalValidLFPs);
figure
h=histogram(wellBurstStartAmp,logspace(1,4,13),"Normalization","probability");
hVals=h.Values;
hBinCenters=convert_edges_2_centers(h.BinEdges);
% xticks(-360:60:360)
xlabel("Axon Theta Amp uV")
ylabel("Soma Burst Start Probability")
ax=gca;
ax.FontSize=16;
xlim([0,2000])
% hold on
% f=fit(hVals',hBinCenters',"log10");
% plot(f,hVals,hBinCenters)
ax.XScale="log";
% ax.YScale="log";

% figure
% plot()
%% Regression tests
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==4 & well_spike_dyn.regi==4);
sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,4,"G10",targetElecs,well_spike_dyn,20,thresh_mult,...
    "D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 24574 210715 21div 210806_1.h5\")