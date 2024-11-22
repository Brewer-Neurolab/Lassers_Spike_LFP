%% Axon to Well Amplitude and Spike Relation Test for Theta

%% Setup
clear
clc
close all

data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\Theta\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
fs=data.re_fs;
data=data.filtered_data;
t_rec=300;
%% plot axon data tagged

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*fs*3; %1 cycle of fastest theta

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2;
minLFPLength=(1/10)*minLFPCycles*fs;

identify_lfps(data,fs,t_rec,minLFPLength,minLFPCycles,nsamples_combine_thresh);