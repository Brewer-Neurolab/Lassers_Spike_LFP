%% Setup
clc;
clear;
close all;

parent_axons_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A";
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

axon_spikes=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
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

all_reg=[interRegions,"EC-CA1"];


%%Draft up before coding
%1. Import axon data
%2. Filter into high/low gamma and theta regions (100-300, 30-100, 4-10)
%3. Hilbert each of these
%4. Find phase of theta and amp of low/high gamma
%5. Find PAC and graph on polarplot

%% Compute theta-gamma regions

for nFF=1:height(ff_axon_tbl)
    % Axon Setup
    data=load(fullfile(parent_axons_dir,axons_folders(ff_axon_tbl.fi(nFF)),ff_axon_tbl.Electrode(nFF)+".mat"));
    data=data.data;
    fs=25000;
    t_rec=300;
    re_fs=1000;
    re_t=0:1/re_fs:t_rec-(1/re_fs);

    %downsample data
    data=resample(data,re_fs,fs);

    %create full sampled time steps
    t=0:1/fs:t_rec-(1/fs);

    %filter for high gamma
    [A,B,C,D]=butter(8,300/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,100/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    high_gamma=filtfilt(sos,g,LFP_filt);

    %filter for low gamma
    [A,B,C,D]=butter(8,100/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,30/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    low_gamma=filtfilt(sos,g,LFP_filt);
    
    %filter for theta
    [A,B,C,D]=butter(8,10/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,4/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    theta=filtfilt(sos,g,LFP_filt);

    HighGammaHilbert=hilbert(high_gamma);
    LowGammaHilbert=hilbert(low_gamma);
    ThetaHilbert=hilbert(theta);

    % %define max number of samples for combining LFPs
    % nsamples_combine_thresh_high=(1/300)*re_fs*3;
    % % nsamples_combine_thresh=[];
    % 
    % %define min lfp length as 2x shortest high gamma cycle
    % minLFPLength_H=(1/300)*minLFPCycles*re_fs;
    % 
    % [LFPEndPts_High,LFPAmplitude_High,LFPHilbert_High]=identify_lfps(high_gamma,re_fs,t_rec,thresh_mult,minLFPLength_H,minLFPCycles,nsamples_combine_thresh_high);
    % LFPAngles=angle(LFPHilbert_High)*(180/pi);
    % 
    % validLFPIndex_H=[];
    % for nEndPts=1:size(LFPEndPts_High,1)
    %     validLFPIndex_H=[validLFPIndex_H,LFPEndPts_High(nEndPts,1):LFPEndPts_High(nEndPts,2)];
    % end
    % validLFPIndex_H=unique(validLFPIndex_H);
    % validLFPIndex_H(validLFPIndex_H<=0 | validLFPIndex_H>length(re_t))=[];
    % logicalValidLFPs_H=zeros(1,length(re_t));
    % logicalValidLFPs_H(validLFPIndex_H(validLFPIndex_H>0 & validLFPIndex_H<length(re_t)))=1;
    % logicalValidLFPs_H=logical(logicalValidLFPs_H);

    thresh_mult=1;
    minLFPCycles=2; %default 2

    %define max number of samples for combining LFPs
    nsamples_combine_thresh_T=(1/10)*re_fs*3;
    minLFPLength_T=(1/10)*minLFPCycles*re_fs;

    [LFPEndPts_T,LFPAmplitude_T,LFPHilbert_T]=identify_lfps(theta,re_fs,t_rec,thresh_mult,minLFPLength_T,minLFPCycles,nsamples_combine_thresh_T);
    % LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));
    LFPAngles=angle(LFPHilbert_T)*(180/pi);

    validLFPIndex_T=[];
    for nEndPts=1:size(LFPEndPts_T,1)
        validLFPIndex_T=[validLFPIndex_T,LFPEndPts_T(nEndPts,1):LFPEndPts_T(nEndPts,2)];
    end
    validLFPIndex_T=unique(validLFPIndex_T);
    validLFPIndex_T(validLFPIndex_T<=0 | validLFPIndex_T>length(re_t))=[];
    logicalValidLFPs_T=zeros(1,length(re_t));
    logicalValidLFPs_T(validLFPIndex_T(validLFPIndex_T>0 & validLFPIndex_T<length(re_t)))=1;
    logicalValidLFPs_T=logical(logicalValidLFPs_T);

    amp_high=abs(HighGammaHilbert);
    amp_high=amp_high(logicalValidLFPs_T);
    amp_low=abs(LowGammaHilbert);
    amp_low=amp_low(logicalValidLFPs_T);
    % amp_t=abs(ThetaHilbert);
    % amp_t=amp_t(logicalValidLFPs_T);
    % phase_h=angle(ThetaHilbert);
    % phase_h=phase_h(logicalValidLFPs_H);
    phase_t=angle(ThetaHilbert);
    phase_t=phase_t(logicalValidLFPs_T);

    PAC_high = abs(mean(amp_high.*exp(1i*phase_t)));
    PAC_low = abs(mean(amp_low.*exp(1i*phase_t)));
    %PAC_theta = abs(mean(amp_t.*exp(1i*phase_t)));
    ff_axon_tbl.PACHigh(nFF) = PAC_high;
    ff_axon_tbl.PACLow(nFF) = PAC_low;
    %%ff_axon_tbl.PACTheta(nFF) = PAC_theta;

    figure
    polarplot(phase_t,amp_high);

    figure
    polarplot(phase_t, amp_low);
    
    disp(nFF)
end


