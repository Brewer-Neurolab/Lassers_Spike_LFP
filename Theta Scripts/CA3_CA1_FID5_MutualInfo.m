clear
clc
close all
%% G12 Axon Setup
data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
data=data.data;
fs=25000;
t_rec=300;
re_fs=1000;
re_t=0:1/re_fs:t_rec-(1/re_fs);

%downsample data
data=resample(data,re_fs,fs);

%create full sampled time steps
t=0:1/fs:t_rec-(1/fs);

%filter for theta
[A,B,C,D]=butter(8,10/(re_fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,4/(re_fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
theta_G12=filtfilt(sos,g,LFP_filt);

%get axon spikes
axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
axon_spikes=axon_spikes.cluster_class(:,2)/1000;
axon_spike_train_12=zeros(1,length(t));
axon_spike_train_G12(ismembertol(t,axon_spikes))=1;

%% G11 Axon Setup
data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\G11.mat");
data=data.data;
fs=25000;
t_rec=300;
re_fs=1000;
re_t=0:1/re_fs:t_rec-(1/re_fs);

%downsample data
data=resample(data,re_fs,fs);

%create full sampled time steps
t=0:1/fs:t_rec-(1/fs);

%filter for theta
[A,B,C,D]=butter(8,10/(re_fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,4/(re_fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
theta_G11=filtfilt(sos,g,LFP_filt);

%get axon spikes
axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G11.mat");
axon_spikes=axon_spikes.cluster_class(:,2)/1000;
axon_spike_train_G11=zeros(1,length(t));
axon_spike_train_G11(ismembertol(t,axon_spikes))=1;

%% Mutual info G12 G11 FID 5

G12Hilbert=hilbert(theta_G12);
G12Amp=abs(G12Hilbert);
G12Angle=angle(G12Hilbert);
G11Hilbert=hilbert(theta_G11);
G11Amp=abs(G11Hilbert);
G11Angle=angle(G11Hilbert);

% Get optimal number of bins TODO

n_bins=20;

figure
edgesG12=logspace(log10(min(G12Amp)),log10(max(G12Amp)),n_bins+1);
hG12=histogram(G12Amp,edgesG12);
probG12=hG12.Values./sum(hG12.Values);
SEG12=-sum(probG12.*log2(probG12+eps));
xscale(gca,"log")
xlim([min(G12Amp),max(G12Amp)])

figure
edgesG11=logspace(log10(min(G11Amp)),log10(max(G11Amp)),n_bins+1);
hG11=histogram(G11Amp,edgesG11);
probG11=hG11.Values./sum(hG11.Values);
SEG11=-sum(probG11.*log2(probG11+eps));
xscale(gca,"log")
xlim([min(G11Amp),max(G11Amp)])

[~,~,bins1]=histcounts(G12Amp,edgesG12);
[~,~,bins2]=histcounts(G11Amp,edgesG11);

jp=zeros(n_bins);
for i1=1:n_bins
    for i2=1:n_bins
        jp(i1,i2)=sum(bins1==i1 & bins2==i2);
    end
end

jp=jp./sum(sum(jp));

JSE=-sum(sum(jp.*log2(jp+eps)));

MutInfo=SEG12+SEG11-JSE;