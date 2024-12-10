%% Axon to Well Amplitude and Spike Relation Test for Theta

%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\Theta\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
re_fs=data.re_fs;
data=data.filtered_data;
t_rec=300;
re_t=0:1/re_fs:t_rec-(1/re_fs);

well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;
%% plot axon data tagged

%define max number of samples for combining LFPs
% nsamples_combine_thresh=(1/10)*re_fs*3; %1 cycle of fastest theta
nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=0.2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,minLFPLength,minLFPCycles,nsamples_combine_thresh);
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

spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\E10_spikes.mat");
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
well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};
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
%% Shuffle control
well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;

MI_Shuffled_BurstStart=[];
for nPerm=1:nPermutations
    % figure
    wellBurstStartAngles=shuffledLFPAngles(nPerm,logicalBurstStarts&logicalValidLFPs);
    phaseBurstStartProb=histogram(wellBurstStartAngles,0:20:360,"Normalization","probability");
    MI_Shuffled_BurstStart(nPerm)=spike_amplitude_MI(phaseBurstStartProb.Values);
    wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
end

figure
h=histogram(MI_Shuffled_BurstStart);
xline(MI_BurstStart,'Color','r')

pValMI_BurstStart=sum(MI_Shuffled_BurstStart>=MI_BurstStart)/length(MI_Shuffled_BurstStart);

xlabel("Shuffled MI")
ylabel("MI Counts")
ax=gca;
ax.FontSize=16;

%% Loop through all target well channels and compare with source axon

MI_Table=table();

CA1_Electrodes=well_spike_dyn.channel_name(well_spike_dyn.fi==6 & well_spike_dyn.regi==4);
row=1;

for nElec=1:length(CA1_Electrodes)
    myElec=CA1_Electrodes(nElec);
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name==myElec};
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
    title(myElec+" "+MI_BurstStart)

    figure
    histogram(wellBurstStartAngles,-360:20:360,"Normalization","probability")
    xticks(-360:60:360)
    xlabel("Axon Theta Angle")
    ylabel("Soma Burst Start Probability")
    ax=gca;
    ax.FontSize=16;
    title(myElec+" MI="+MI_BurstStart)

    % check for significance
    MI_Shuffled_BurstStart=[];
    for nPerm=1:nPermutations
        % figure
        wellBurstStartAngles=shuffledLFPAngles(nPerm,logicalBurstStarts&logicalValidLFPs);
        phaseBurstStartProb=histogram(wellBurstStartAngles,0:20:360,"Normalization","probability");
        MI_Shuffled_BurstStart(nPerm)=spike_amplitude_MI(phaseBurstStartProb.Values);
        wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    end

    figure
    h=histogram(MI_Shuffled_BurstStart);
    xline(MI_BurstStart,'Color','r')

    pValMI_BurstStart=sum(MI_Shuffled_BurstStart>=MI_BurstStart)/length(MI_Shuffled_BurstStart);

    xlabel("Shuffled MI")
    ylabel("MI Counts")
    ax=gca;
    ax.FontSize=16;

    title(myElec+" pval="+pValMI_BurstStart)

    %update MI table
    MI_Table.fi(row)=6;
    MI_Table.Tunnel_Electrode(row)="G12";
    MI_Table.Tunnel_reg(row)="CA3-CA1";
    MI_Table.is_ff(row)=1;
    MI_Table.Well_Electrode(row)=CA1_Electrodes(nElec);
    MI_Table.Well_reg(row)="CA1";
    MI_Table.MI(row)=MI_BurstStart;
    MI_Table.pval(row)=pValMI_BurstStart;

    row=row+1;
end

%% Well spikes in bursts angles at high axon LFP
well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};
well_burst_bounds=round(remap(well_burst_bounds,1,length(t),1,length(re_t)));
logicalIsBursting=zeros(1,length(re_t));

for nBurst=1:size(well_burst_bounds,1)
    logicalIsBursting(well_burst_bounds(nBurst,1):well_burst_bounds(nBurst,2))=1;
end

wellBurstSpikeAngles=LFPAngles(logicalIsBursting&logicalValidLFPs);
wellBurstSpikeAngles=[wellBurstSpikeAngles-360,wellBurstSpikeAngles];

figure
histogram(wellBurstSpikeAngles,-360:20:360)
xticks(-360:60:360)
xlabel("Axon Theta Angle")
ylabel("Soma Burst Spike Count")
ax=gca;
ax.FontSize=16;

figure
histogram(wellBurstSpikeAngles,-360:20:360,"Normalization","probability")
xticks(-360:60:360)
xlabel("Axon Theta Angle")
ylabel("Soma Burst Spike Probability")
ax=gca;
ax.FontSize=16;
%% compare spikes at LFP

wellSpikeAmp=LFPAmplitude(logicalSpikes&logicalValidLFPs);

logbinAmp=logspace(1,4,13);

figure
histogram(wellSpikeAmp,logbinAmp)
xlabel("Axon Theta Amplitude")
ylabel("Soma Spike")

ax=gca;
ax.XScale="log";

%% Spikes in bursts in LFP after burst start

%find burst starts
well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;

%define burst ends
well_burst_ends=well_burst_bounds(:,2);
% remap burst starts to new sampling rate
well_burst_ends=round(remap(well_burst_ends,1,length(t),1,length(re_t)));
logicalBurstEnds=zeros(1,length(re_t));
logicalBurstEnds(well_burst_ends)=1;

%find burst starts/ends in high LFP
highAmpBurstStarts=logicalBurstStarts(logicalBurstStarts & logicalValidLFPs);

%create new bursting start/end matrix
burstIdx=round(well_burst_starts)==round(find(highAmpBurstStarts));
highBursts=well_burst_starts(burstIdx,:);

spikesInHighBurstsAngles=[];
for nBursts=1:size(highBursts,1)
    validBurstIdx=highBursts(nBursts,1):highBursts(nBursts,2);
end
%% 3D Graph of E10

well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};
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

X=[wellBurstStartAngles;wellBurstStartAmp]';
edges={[0:40:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10)};

hist3(X,'Edges',edges)
xlabel("Axon Phase Angle")
ylabel("Axon Amplitude")
zlabel("Soma Burst Starts Counts")

ax=gca;
ax.YScale="log";