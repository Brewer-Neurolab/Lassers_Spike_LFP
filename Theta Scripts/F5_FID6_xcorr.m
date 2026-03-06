%% Xcorr EC-CA3 FID 6 F5
clear
clc
close all
%load downsampled LFP
LFP=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\4x 33168 210715 21div 210806_1_mat_files\F5.mat");
LFP=LFP.re_LFP;
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\times_F5.mat");
spikes=spikes.cluster_class(:,2);
%5SD spikes
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

% EC_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==1 & well_spike_dyn.fi==6);
EC_Elec=["A4","A5","B3","B4","B5","C2","C3","C4","C5","D1","D2","D3","D4","D5","E1","E2","E3","E4","E5"]

fs=25000;
re_fs=1000;
t=1/fs:1/fs:300;
re_t=1/re_fs:1/re_fs:300;
t_rec=300;
%% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(LFP,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
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
logicalValidLFPs=logical(logicalValidLFPs);


%% xcorr CA3 LFP
xcorr_table=table();

for nElec=1:length(EC_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled wells\Theta\4x 33168 210715 21div 210806_1.h5\"+EC_Elec(nElec)+".mat");
    wellLFP=myElec.filtered_data;

    [r,l]=xcorr(LFP,wellLFP,re_fs*1,"normalized");
    % [r,l]=xcorr(LFP(logicalValidLFPs),wellLFP(logicalValidLFPs),re_fs*2,"normalized");
    figure
    plot(l/re_fs*1000,r,"LineWidth",2)
    title(EC_Elec(nElec))

    [~,P]=corrcoef(wellLFP,LFP);
    % [~,P]=corrcoef(wellLFP(logicalValidLFPs),LFP(logicalValidLFPs));

    xcorr_table.elec(nElec)=EC_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.rMax{nElec}=max(abs(r));
    xcorr_table.l{nElec}=l;
    xcorr_table.lMax{nElec}=l(abs(r)==max(abs(r)));
    xcorr_table.pval(nElec)=P(1,2);
    xlabel("Lag (ms)")
    ylabel("r normalized")
    set(gca,"FontSize",40)
    xlim([-200,200])

    ax=gca;
    ax.LineWidth=4;
    ax.TickLength=[0.05 0.05];
    axis square
end