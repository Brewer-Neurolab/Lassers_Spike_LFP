clear
clc

%% Axon Setup
LFP=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
LFP=LFP.data;
fs=25000;
t_rec=300;

t=0:1/fs:t_rec-(1/fs);

[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,LFP);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,LFP_filt);

%for derivative of LFP data
% LFP_filt=[diff(LFP_filt),0];

% G12_env=abs(hilbert(LFP_filt));
%% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(LFP_filt,fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

LFP_filt(~logicalValidLFPs)=0;
LFPAmplitude(~logicalValidLFPs)=0;
LFPAngles(~logicalValidLFPs)=0;
%% E12 setup

soma=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10.mat");
soma=soma.data;

spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;

[b,a]=ellip(2,0.1,40,[300,3000]/(fs/2));
filtSpikeData=filtfilt(b,a,soma);

% spikeHilbert=abs(hilbert(filtSpikeData));
spikeHilbert=envelope(abs(hilbert(filtSpikeData)),0.010*fs,'peak'); % create a window of 150 ms for minimum ISI in a burst to capture bursts rather than individual spikes

% spikeHilbert=abs(hilbert(spike_train));
% spikeHilbert=envelope(spike_train,5*fs,'analytic');
%% plot both
figure
plot(t,LFP_filt)
hold on
plot(t,LFPAmplitude)
hold off

figure
plot(t,filtSpikeData)
% plot(t,spike_train)
hold on
plot(t,spikeHilbert)
hold off

figure
tiledlayout(2,1)
nexttile
plot(t,LFP_filt)
hold on
plot(t,LFPAmplitude)
hold off
ax1=gca;

nexttile
plot(t,filtSpikeData)
% plot(t,spike_train)
hold on
plot(t,spikeHilbert)
hold off
ax2=gca;

linkaxes([ax1,ax2],'x')
%% xcorr angles

[r,l]=xcorr(spikeHilbert,LFPAngles,1*fs,"normalized");
figure
plot(l/fs,r)
xlabel("Seconds")
ylabel("r")
ax=gca;
ax.FontSize=18;

%% xcorr ampsXpeaks

[r,l]=xcorr(spikeHilbert,LFPAmplitude,1*fs,"normalized");
figure
plot(l/fs,r)
xlabel("Seconds")
ylabel("r")
ax=gca;
ax.FontSize=18;

%% xcorr with all well spikes
parent_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\";
%3.5SD spikes
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

CA1_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==4 & well_spike_dyn.fi==5);

for nElec=1:4%length(CA1_Elec)
    soma=load(parent_dir+CA1_Elec(nElec));
    soma=soma.data;

    [b,a]=ellip(2,0.1,40,[300,3000]/(fs/2));
    filtSpikeData=filtfilt(b,a,soma);

    % spikeHilbert=abs(hilbert(filtSpikeData));
    spikeHilbert=envelope(abs(hilbert(filtSpikeData)),0.010*fs,'peak'); % create a window of 150 ms for minimum ISI in a burst to capture bursts rather than individual spikes

    figure('Name',CA1_Elec(nElec),'NumberTitle','off')
    tFig=tiledlayout(2,1);
    nexttile
    plot(t,LFP_filt,'k')
    hold on
    plot(t,LFPAmplitude,'g')
    ylabel("uV")
    yyaxis right
    plot(t,LFPAngles,'m')
    ylabel("Degrees")
    hold off
    xlabel("Seconds")
    ax1=gca;
    set(gca,'FontSize',18)
    ax1.YAxis(1).Color='k';
    ax1.YAxis(2).Color='m';
    
    nexttile
    plot(t,filtSpikeData,'k')
    ylabel("uV")
    % plot(t,spike_train)
    hold on
    plot(t,spikeHilbert,'g')
    hold off
    xlabel("Seconds")
    ax2=gca;
    set(gca,'FontSize',18)

    linkaxes([ax1,ax2],'x')

    % title(tFig,CA1_Elec(nElec))

    [r,l]=xcorr(spikeHilbert,LFPAngles,1*fs,"normalized");
    figure('Name',CA1_Elec(nElec)+" XCorr",'NumberTitle','off')
    plot(l/fs,r)
    xlabel("Seconds")
    ylabel("r")
    title(CA1_Elec(nElec))
    ax=gca;
    ax.FontSize=18;

    [r,l]=xcorr(spikeHilbert,wrapTo360(LFPAngles+180),1*fs,"normalized");
    figure('Name',CA1_Elec(nElec)+" XCorr+180",'NumberTitle','off')
    plot(l/fs,r)
    xlabel("Seconds")
    ylabel("r")
    title(CA1_Elec(nElec))
    ax=gca;
    ax.FontSize=18;
    
end

