%% Axon to Well Amplitude and Spike Relation Test for Theta

%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\Theta\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
re_fs=data.re_fs;
data=data.filtered_data;
t_rec=300;
re_t=0:1/re_fs:t_rec-(1/re_fs);
fs=25000;
t=0:1/fs:t_rec-(1/fs);

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
shuffledLFPAngles=[];
for nPerm=1:100
    randPermAngles=randperm(length(LFPAngles));
    shuffledLFPAngles(nPerm,:)=LFPAngles(randPermAngles);
end

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
end
logicalValidLFPs=zeros(1,length(re_t));
logicalValidLFPs(validLFPIndex)=1;

%% Loop through all target well channels and compare with source axon

MI_Table=table();

CA3_Electrodes=well_spike_dyn.channel_name(well_spike_dyn.fi==5 & well_spike_dyn.regi==3);
row=1;

for nElec=1:length(CA3_Electrodes)
    myElec=CA3_Electrodes(nElec);
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
    for nPerm=1:100
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
    MI_Table.Well_Electrode(row)=CA3_Electrodes(nElec);
    MI_Table.Well_reg(row)="CA1";
    MI_Table.MI(row)=MI_BurstStart;
    MI_Table.pval(row)=pValMI_BurstStart;

    row=row+1;
end




