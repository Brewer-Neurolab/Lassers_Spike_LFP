%% Axon LFP->Axon Spike Amplitude correlations

clear
clc
close all
%% Load and filter axon LFP for spikes
% data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
fs=data.sr;
data=data.data;

% [A,B,C,D]=butter(8,3000/(fs/2),'low');
% [sos,g]=ss2sos(A,B,C,D);
% spike_data=filtfilt(sos,g,data);
% 
% [A,B,C,D]=butter(8,300/(fs/2),'high');
% [sos,g]=ss2sos(A,B,C,D);
% spike_data=filtfilt(sos,g,spike_data);

[b,a]=ellip(2,0.1,40,[300,3000]/(fs/2));
filtSpikeData=filtfilt(b,a,data);

[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtLFPdata=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtLFPdata=filtfilt(sos,g,filtLFPdata);

LFPAmp=abs(hilbert(filtLFPdata));

% t=1/fs:1/fs:300;
t=0:1/fs:300-(1/fs);

figure
plot(t,filtLFPdata)

%% Get Spike Amplitudes
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\times_G12.mat"); 
spikes=spikes.cluster_class(:,2)/1000; %ms to s
[spike_idx,spikes_used_idx]=ismembertol(t,spikes);

spike_amps=filtSpikeData(spike_idx);

figure
plot(t,filtSpikeData)
hold on
plot(t(spike_idx),filtSpikeData(spike_idx),'r.')
hold off

%% Loop through wells and get pearson and linear regression
% parent_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\";
% active_electrodes
log_bin_amps=logspace(1,log10(max(LFPAmp(spike_idx))),15);
figure
scatter(LFPAmp(spike_idx),spike_amps)
xlabel("Axon Theta LFP uV")
ylabel("Spike Amplitude uV")

x=LFPAmp(spike_idx);
x=x(spike_amps<-800);
y=spike_amps(spike_amps<-800);
figure
scatter(x,y)
xlabel("Axon Theta LFP uV")
ylabel("Spike Amplitude uV")
ax=gca;
ax.XScale="log";
% ax.YScale="log";
set(gca,"FontSize",18)
ax.TickLength(1)=0.05;
xlim([min(log_bin_amps),max(log_bin_amps)])

%% Histogram of Spiking
% log_bin_amps=logspace(log10(min(LFPAmp(spike_idx))),log10(max(LFPAmp(spike_idx))),15);

figure
histogram(LFPAmp(spike_idx),log_bin_amps)
ax=gca;
ax.XScale="log";
xlabel("Axon Theta Amp uV")
ylabel("Spike counts")
set(gca,"FontSize",18)
ax.TickLength(1)=0.05;