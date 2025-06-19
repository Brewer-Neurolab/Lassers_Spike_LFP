%% Xcorr CA3-CA1 FID 5 G11
clear
clc
close all
%load downsampled LFP
LFP=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\Theta\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
LFP=LFP.filtered_data;
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
spikes=spikes.cluster_class(:,2);
%5SD spikes
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

CA1_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==4 & well_spike_dyn.fi==5);
CA3_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==3 & well_spike_dyn.fi==5);

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
%% xcorr CA1 Spikes
xcorr_table=table();

% try for 50 bins on x axis, downsample

spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes/1000))=1;

re_spike_train=ceil(remap(find(spike_train),1,length(t),1,length(re_t)));
myZeros=zeros(1,length(re_t));
myZeros(re_spike_train)=1;
spike_train=myZeros;

for nElec=1:length(CA1_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\"+CA1_Elec(nElec)+"_spikes.mat");
    my_spike_train=zeros(1,length(t));
    my_spike_train(ismembertol(t,myElec.index/1000))=1;

    re_my_spike_train=ceil(remap(find(my_spike_train),1,length(t),1,length(re_t)));
    myZeros=zeros(1,length(re_t));
    myZeros(re_my_spike_train)=1;
    my_spike_train=myZeros;

    [r,l]=xcorr(spike_train,my_spike_train,re_fs*0.01,"normalized");
    figure
    plot(l/re_fs*1000,r)
    title(CA1_Elec(nElec))

    [~,P]=corrcoef(spike_train,my_spike_train);

    xcorr_table.elec(nElec)=CA1_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.l{nElec}=l;
    xcorr_table.pval(nElec)=P(1,2);
end
%% xcorr CA3 Spikes
xcorr_table=table();

% try for 50 bins on x axis, downsample

spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes/1000))=1;

re_spike_train=ceil(remap(find(spike_train),1,length(t),1,length(re_t)));
myZeros=zeros(1,length(re_t));
myZeros(re_spike_train)=1;
spike_train=myZeros;

for nElec=1:length(CA3_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\"+CA3_Elec(nElec)+"_spikes.mat");
    my_spike_train=zeros(1,length(t));
    my_spike_train(ismembertol(t,myElec.index/1000))=1;

    re_my_spike_train=ceil(remap(find(my_spike_train),1,length(t),1,length(re_t)));
    myZeros=zeros(1,length(re_t));
    myZeros(re_my_spike_train)=1;
    my_spike_train=myZeros;

    [r,l]=xcorr(spike_train,my_spike_train,re_fs*0.01,"normalized");
    figure
    plot(l/re_fs*1000,r)
    title(CA3_Elec(nElec))

    [~,P]=corrcoef(spike_train,my_spike_train);

    xcorr_table.elec(nElec)=CA3_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.l{nElec}=l;
    xcorr_table.pval(nElec)=P(1,2);
end
%% xcorr CA1 LFP
xcorr_table=table();

for nElec=1:length(CA1_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled wells\Theta\4x 33152 210715 21div 210806_1.h5\"+CA1_Elec(nElec)+".mat");
    wellLFP=myElec.re_LFP;

    % [r,l]=xcorr(wellLFP,LFP,re_fs*0.1,"normalized");
    [r,l]=xcorr(wellLFP(logicalValidLFPs),LFP(logicalValidLFPs),re_fs*0.1,"normalized");
    figure
    plot(l/re_fs*1000,r)
    title(CA1_Elec(nElec))

    % [~,P]=corrcoef(LFP,wellLFP);
    [~,P]=corrcoef(LFP(logicalValidLFPs),wellLFP(logicalValidLFPs));

    xcorr_table.elec(nElec)=CA1_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.rMax{nElec}=max(abs(r));
    xcorr_table.l{nElec}=l;
    xcorr_table.lMax{nElec}=l(abs(r)==max(abs(r)));
    xcorr_table.pval(nElec)=P(1,2);
    xlabel("Lag (ms)")
    ylabel("r normalized")
    set(gca,"FontSize",18)
end
%% xcorr CA3 LFP
xcorr_table=table();

for nElec=1:length(CA3_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled wells\Theta\4x 33152 210715 21div 210806_1.h5\"+CA3_Elec(nElec)+".mat");
    wellLFP=myElec.re_LFP;

    [r,l]=xcorr(LFP,wellLFP,re_fs*1,"normalized");
    % [r,l]=xcorr(LFP(logicalValidLFPs),wellLFP(logicalValidLFPs),re_fs*2,"normalized");
    figure
    plot(l/re_fs*1000,r,"LineWidth",2)
    title(CA3_Elec(nElec))

    [~,P]=corrcoef(wellLFP,LFP);
    % [~,P]=corrcoef(wellLFP(logicalValidLFPs),LFP(logicalValidLFPs));

    xcorr_table.elec(nElec)=CA3_Elec(nElec);
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