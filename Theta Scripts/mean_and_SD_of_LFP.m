%% Axon to Well Amplitude and Spike Relation Test for Theta
%% Setup
clear
clc
close all

saveDir="C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Theta Scripts";

parent_axons_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\theta";
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

parent_wells_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min";
wells_dir=dir(parent_wells_dir);
wells_folders=string({wells_dir.name});
wells_folders=wells_folders([wells_dir.isdir]);
wells_folders=wells_folders(3:end);

axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axon_spikes=axon_spikes.allregion_unit_matched_stim;

interRegions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];
subregions=["EC","DG","CA3","CA1"];
ff_axon_tbl=table();

row=1;
for fi=1:length(axon_spikes)
    for nelec=1:height(axon_spikes{fi})
        if ~isempty(axon_spikes{fi}.up_ff{nelec}) | (~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC")
            ff_axon_tbl.fi(row)=fi;
            ff_axon_tbl.Subregion(row)=axon_spikes{fi}.Subregion(nelec);
            ff_axon_tbl.interRegi(row)=find(interRegions==axon_spikes{fi}.Subregion(nelec));
            electrodes=axon_spikes{fi}.("Electrode Pairs")(nelec);
            electrodes=split(electrodes,{'-'});
            ff_axon_tbl.Electrode(row)=electrodes(1);
            tunnelReg=split(axon_spikes{fi}.Subregion(nelec),{'-'});
            ff_axon_tbl.FFReg(row)=tunnelReg(2);
            ff_axon_tbl.subi(row)=find(subregions==tunnelReg(2));

            if isempty(axon_spikes{fi}.up_ff{nelec}) & ~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC"
                ff_axon_tbl.Subregion(row)="EC-CA1";
                ff_axon_tbl.FFReg(row)="CA1";
                ff_axon_tbl.subi(row)=4;
            end

            row=row+1;
        end
    end
end

%3.5 SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

all_reg=[interRegions,"EC-CA1"];
%% Compute Mean and SD tables
clc

%testing range
testing_idx=[find(ff_axon_tbl.Subregion=="CA3-CA1")]';

lfpPropsTab=ff_axon_tbl;

for nFF=1:height(ff_axon_tbl)
    data=load(fullfile(parent_axons_dir,axons_folders(ff_axon_tbl.fi(nFF)),ff_axon_tbl.Electrode(nFF)+".mat"));
    % re_fs=data.re_fs;
    re_fs=1000;
    fs=25000;
    data=data.filtered_data;
    t_rec=300;
    re_t=0:1/re_fs:t_rec-(1/re_fs);
    t=0:1/fs:t_rec-(1/fs);

    % plot axon data tagged
    thresh_mult=1;

    %define max number of samples for combining LFPs
    nsamples_combine_thresh=(1/10)*re_fs*3;
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest theta cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/10)*minLFPCycles*re_fs;

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    validLFPIndex=[];
    for nEndPts=1:size(LFPEndPts,1)
        validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
    end
    logicalValidLFPs=zeros(1,length(re_t)); % uncomment for lower bound LFP
    % logicalValidLFPs=ones(1,length(re_t)); % considers all LFPs
    logicalValidLFPs(validLFPIndex)=1;

    % propsTabVars=ff_axon_tbl.Variables
    % lfpPropsTab(nFF,1:6)=ff_axon_tbl(nFF,:);
    lfpPropsTab.RMS(nFF)=rms(data);
    lfpPropsTab.MeanLFPAmp_uV(nFF)=mean(LFPAmplitude);
    lfpPropsTab.LFPAmp_SD(nFF)=std(LFPAmplitude);
    lfpPropsTab.("perc_above_mean+sd")(nFF)=(sum(logicalValidLFPs)/length(logicalValidLFPs))*100;
    lfpPropsTab.perc_above_5uV(nFF)=(sum(LFPAmplitude>=5)/length(LFPAmplitude))*100;
end

writetable(lfpPropsTab,fullfile("C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Theta Scripts","lfpPropsTab.xlsx"))

%% Plot by subregion

% subreg_mean=[];
% subreg_meanSD=[];
% subreg_SD=[];
lfp_stats_tbl=table();

for i=1:length(all_reg)
    lfp_stats_tbl.subregion(i)=all_reg(i);
    lfp_stats_tbl.subreg_mean(i)=mean(lfpPropsTab.MeanLFPAmp_uV(lfpPropsTab.Subregion==all_reg(i)));
    lfp_stats_tbl.subreg_meanSD(i)=mean(lfpPropsTab.LFPAmp_SD(lfpPropsTab.Subregion==all_reg(i)));
    lfp_stats_tbl.subreg_SD(i)=std(lfpPropsTab.MeanLFPAmp_uV(lfpPropsTab.Subregion==all_reg(i)));
    lfp_stats_tbl.mean_rms(i)=mean(lfpPropsTab.RMS(lfpPropsTab.Subregion==all_reg(i)));
end

writetable(lfp_stats_tbl,fullfile("C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Theta Scripts","lfpStatsTab.xlsx"))

%% Has target?

