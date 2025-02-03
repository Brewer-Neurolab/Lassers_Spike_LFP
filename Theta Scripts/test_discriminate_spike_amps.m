%% FID 5 M9 Spike amplitude histogram
clear
clc
fs=25000;
t=1/fs:1/fs:300;

M9_signal=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\M9.mat");
M9_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\M9_spikes.mat");

d=designfilt("bandpassiir",FilterOrder=6,SampleRate=fs,HalfPowerFrequency1=300,HalfPowerFrequency2=3000);

M9_signal=filtfilt(d,M9_signal.data);

M9_spike_train=zeros(1,length(t));

M9_spike_train(ismembertol(t,M9_spikes.index/1000))=1;

figure
histogram(M9_signal(logical(M9_spike_train)))