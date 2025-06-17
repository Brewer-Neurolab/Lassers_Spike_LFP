%% Axon to Well Amplitude and Spike Relation Test for Theta

%% Setup
clear
clc
close all

data=load("C:\Users\ssk78\Desktop\Brewer LFP\Tunnels\High_Gamma\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
re_fs=data.re_fs;
=======
%re_fs=data.re_fs;
re_fs = 1000;
>>>>>>> Stashed changes
=======
%re_fs=data.re_fs;
re_fs = 1000;
>>>>>>> Stashed changes
=======
%re_fs=data.re_fs;
re_fs = 1000;
>>>>>>> Stashed changes
data=data.filtered_data;


t_rec=300;
re_t=0:1/re_fs:t_rec-(1/re_fs);

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
well_spike_dyn=load("C:\Users\ssk78\Desktop\Brewer LFP\well_spike_dynamics_table_hfs_3-5.mat");
=======
well_spike_dyn=load("C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
>>>>>>> Stashed changes
=======
well_spike_dyn=load("C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
>>>>>>> Stashed changes
=======
well_spike_dyn=load("C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
>>>>>>> Stashed changes
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;
%% plot axon data tagged
thresh = 1;
%define max number of samples for combining LFPs
% nsamples_combine_thresh=(1/10)*re_fs*3; %1 cycle of fastest theta
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
nsamples_combine_thresh=0.01*re_fs;
=======
nsamples_combine_thresh=[re_fs*0.05];
>>>>>>> Stashed changes
=======
nsamples_combine_thresh=[re_fs*0.05];
>>>>>>> Stashed changes
=======
nsamples_combine_thresh=[re_fs*0.05];
>>>>>>> Stashed changes

%define min lfp length as 2x shortest theta cycle
minLFPCycles=3; %default 2
minLFPLength=(1/300)*minLFPCycles*re_fs;

<<<<<<< Updated upstream
[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,0.75, minLFPLength,minLFPCycles,nsamples_combine_thresh);
=======
[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh,minLFPLength,minLFPCycles,nsamples_combine_thresh);
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
spikes=load("C:\Users\ssk78\Desktop\Brewer LFP\E10_spikes.mat");
=======
spikes=load("C:\Users\ssk78\Downloads\E10_spikes.mat");
>>>>>>> Stashed changes
=======
spikes=load("C:\Users\ssk78\Downloads\E10_spikes.mat");
>>>>>>> Stashed changes
=======
spikes=load("C:\Users\ssk78\Downloads\E10_spikes.mat");
>>>>>>> Stashed changes
spikes=spikes.index;
fs=25000;
t=0:1/fs:t_rec-(1/fs);


spikes=round(remap(spikes,1,length(t),1,length(re_t)));
logicalSpikes=zeros(1,length(re_t));
logicalSpikes(spikes)=1;

wellSpikeAngles=LFPAngles(logicalSpikes&logicalValidLFPs);
wellSpikeAngles=[wellSpikeAngles-360,wellSpikeAngles];

%% Cummulative Hig Amp Burst start
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==6 & well_spike_dyn.regi==4);
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 9, 0.75)

%% Cummulative Low Amp Burst start
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==6 & well_spike_dyn.regi==4);
cummulative_axon_well_burst_start(t,re_t,~logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn)
=======
sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 20, thresh,"C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\")
>>>>>>> Stashed changes
=======
cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 20, thresh)
%% Cummulative Low Amp Burst start
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==6 & well_spike_dyn.regi==4);
cummulative_axon_well_burst_start(t,re_t,~logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 20, thresh,"C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\")
>>>>>>> Stashed changes
=======
cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 20, thresh)
%% Cummulative Low Amp Burst start
targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==6 & well_spike_dyn.regi==4);
cummulative_axon_well_burst_start(t,re_t,~logicalValidLFPs,LFPAmplitude,LFPAngles,6,"G12",targetElecs,well_spike_dyn, 20, thresh,"C:\Users\ssk78\UC Irvine\BME-Brewer Lab-SPO - Documents\Sam Lassers\Data\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\")
>>>>>>> Stashed changes
