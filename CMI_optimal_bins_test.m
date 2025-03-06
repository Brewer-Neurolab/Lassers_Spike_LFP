% Testing for optimal number of bins

%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\theta\4x 33152 210715 21div 210806_1_mat_files\G10.mat");
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

%define max number of samples for combining LFPs
% nsamples_combine_thresh=(1/10)*re_fs*3; %1 cycle of fastest theta
nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=0.2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,1,minLFPLength,minLFPCycles,nsamples_combine_thresh);
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

%% 3D Graph of E10 AT START OF WELL BURST ONLY PHASE AND AMP SPIKES PER BURST CUMMULATIVE DIST
%cummulative dist may be most accurate

well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==5 & well_spike_dyn.channel_name=="M9"};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;

figure
wellBurstStartAngles=LFPAngles(logicalBurstStarts & logicalValidLFPs);
% wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & logicalValidLFPs);

thetaAmpThresh=std(LFPAmplitude);

%repeat for spikes per burst
repwellBurstStartAngles=[];
repwellBurstStartAmp=[];
repwellSPB=[];

burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
burstIdx=find(burstIdx);
well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==5 & well_spike_dyn.channel_name=="M9"};

for nBursts=1:length(wellBurstStartAngles)
    repwellBurstStartAngles=[repwellBurstStartAngles,wellBurstStartAngles(nBursts)];
    repwellBurstStartAmp=[repwellBurstStartAmp,wellBurstStartAmp(nBursts)];
    repwellSPB=[repwellSPB,well_spb(burstIdx(nBursts))];
end

X=[repwellBurstStartAngles;repwellBurstStartAmp]';
edges={[0:40:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10)};

hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp')
xlabel("Axon Phase Angle (deg)")
ylabel("Axon Burst Start Amplitude (uV)")
zlabel("Cummulative Soma Spikes Per Burst")

xlim([0,360])
yticks(0:10:100)
ylim([(thetaAmpThresh),max(LFPAmplitude)])

ax=gca;
ax.YScale="log";

%% Determine bins

nAngleBins=freedmandiaconis(repwellBurstStartAngles);
nAmpBins=freedmandiaconis(repwellBurstStartAmp);
nSPBBins=freedmandiaconis(log10(repwellSPB));