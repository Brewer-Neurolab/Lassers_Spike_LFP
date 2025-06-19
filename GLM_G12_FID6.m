clear
clc

%% Axon Setup
data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
data=data.data;
fs=25000;
t_rec=300;
re_fs=200;
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
theta=filtfilt(sos,g,LFP_filt);

%filter for delta
[A,B,C,D]=butter(8,4/(re_fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,0.5/(re_fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
delta=filtfilt(sos,g,LFP_filt);

%get axon spikes
axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\times_G12.mat");
axon_spikes=axon_spikes.cluster_class(:,2)/1000;
axon_spike_train=zeros(1,length(t));
axon_spike_train(ismembertol(t,axon_spikes))=1;

%% E12 setup

soma=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\E12.mat");
soma=soma.data;

spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\E12_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;

[b,a]=ellip(2,0.1,40,[300,3000]/(fs/2));
filtSpikeData=filtfilt(b,a,soma);

% spikeHilbert=abs(hilbert(filtSpikeData));
spikeHilbert=envelope(abs(hilbert(filtSpikeData)),0.010*fs,'peak'); % create a window of 150 ms for minimum ISI in a burst to capture bursts rather than individual spikes

% spikeHilbert=abs(hilbert(spike_train));
% spikeHilbert=envelope(spike_train,5*fs,'analytic');

%% Create GLM tbl

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAngle=angle(ThetaHilbert);

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAngle=angle(DeltaHilbert);

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;

glmTbl.WellSpikes=logical(spike_train');
glmTbl.ThetaAmp=ThetaAmp';
glmTbl.ThetaAngle=ThetaAngle';
% glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.DeltaAngle=DeltaAngle';
axon_spike_idx=round(remap(find(axon_spike_train),1,length(t),1,300*re_fs));
axon_spike_train=zeros(1,300*re_fs);
axon_spike_train(axon_spike_idx)=1;
glmTbl.AxonSpikes=logical(axon_spike_train');

modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + ' ...
    'DeltaAmp + DeltaAngle+ AxonSpikes + '...
    'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*AxonSpikes'];
mdl=fitglm(glmTbl,modelspec,'Distribution','binomial');

%% GLM on theta regions zeroed
% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAmp(~logicalValidLFPs)=0;
ThetaAngle=angle(ThetaHilbert);
ThetaAngle(~logicalValidLFPs)=0;

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAmp(~logicalValidLFPs)=0;
DeltaAngle=angle(DeltaHilbert);
DeltaAngle(~logicalValidLFPs)=0;

% glmTbl=table();
% glmTbl.WellSpikes=spike_train';
% glmTbl.ThetaAmp=ThetaAmp';
% glmTbl.ThetaAngle=ThetaAngle';
% glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=DeltaAngle';
% glmTbl.AxonSpikes=axon_spike_train';

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;
glmTbl.WellSpikes=logical(spike_train');
glmTbl.ThetaAmp=resample(ThetaAmp',re_fs,fs);
glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges);
glmTbl.DeltaAmp=resample(DeltaAmp',re_fs,fs);
glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges);
axon_spike_idx=round(remap(find(axon_spike_train),1,length(t),1,300*re_fs));
axon_spike_train=zeros(1,300*re_fs);
axon_spike_train(axon_spike_idx)=1;
glmTbl.AxonSpikes=logical(axon_spike_train');

modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + ' ...
    'DeltaAmp + DeltaAngle+ AxonSpikes + '...
    'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*AxonSpikes'];
mdl_zeroed=fitglm(glmTbl,modelspec,'Distribution','binomial');

%% GLM on theta regions deletion
% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAmp(~logicalValidLFPs)=[];
ThetaAngle=angle(ThetaHilbert);
ThetaAngle(~logicalValidLFPs)=[];

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAmp(~logicalValidLFPs)=[];
DeltaAngle=angle(DeltaHilbert);
DeltaAngle(~logicalValidLFPs)=[];

% glmTbl=table();
% spike_train_del=spike_train;
% spike_train_del(~logicalValidLFPs)=[];
% glmTbl.WellSpikes=spike_train_del';
% glmTbl.ThetaAmp=ThetaAmp';
% glmTbl.ThetaAngle=ThetaAngle';
% glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=DeltaAngle';
% axon_del=axon_spike_train;
% axon_del(~logicalValidLFPs)=[];
% glmTbl.AxonSpikes=axon_del';

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;
glmTbl.WellSpikes=spike_train';
glmTbl.ThetaAmp=resample(ThetaAmp',re_fs,fs);
glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges);
glmTbl.DeltaAmp=resample(DeltaAmp',re_fs,fs);
glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges);
axon_spike_idx=round(remap(find(axon_spike_train),1,length(t),1,300*re_fs));
axon_spike_train=zeros(1,300*re_fs);
axon_spike_train(axon_spike_idx)=1;
glmTbl.AxonSpikes=axon_spike_train';

modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + ' ...
    'DeltaAmp + DeltaAngle+ AxonSpikes + '...
    'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*AxonSpikes'];
mdl_del=fitglm(glmTbl,modelspec,'Distribution','binomial');

%% GLM AxonLFP-Well Spikes NaN outside theta regions

% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(re_t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAmp(~logicalValidLFPs)=NaN;
ThetaAngle=angle(ThetaHilbert);
ThetaAngle(~logicalValidLFPs)=NaN;

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAmp(~logicalValidLFPs)=NaN;
DeltaAngle=angle(DeltaHilbert);
DeltaAngle(~logicalValidLFPs)=NaN;

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;
glmTbl.WellSpikes=logical(spike_train');
glmTbl.ThetaAmp=ThetaAmp';
% glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.ThetaAngle=ThetaAngle';
glmTbl.SinThetaAngle=sin(glmTbl.ThetaAngle);
glmTbl.CosThetaAngle=cos(glmTbl.ThetaAngle);
glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.DeltaAngle=DeltaAngle';
glmTbl.SinDeltaAngle=sin(glmTbl.DeltaAngle);
glmTbl.CosDeltaAngle=cos(glmTbl.DeltaAngle);

% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'SinThetaAngle + CosThetaAngle + SinDeltaAngle + CosDeltaAngle+'...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*'...
%     'SinThetaAngle*CosThetaAngle*SinDeltaAngle*CosDeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + '...
%     'ThetaAmp*ThetaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
modelspec='WellSpikes ~ ThetaAmp*(ThetaAngle+SinThetaAngle+CosThetaAngle)';
% modelspec=['WellSpikes ~ ThetaAmp*ThetaAngle*SinThetaAngle*CosThetaAngle' ...
%     '*DeltaAmp*DeltaAngle*SinDeltaAngle*CosDeltaAngle'];

mdl_zeroed=fitglm(glmTbl,modelspec,'Distribution','binomial','Link','logit');

%% GLM AxonLFP-Well Spikes NaN outside theta regions Moving Sum

% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(re_t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAmp(~logicalValidLFPs)=NaN;
ThetaAngle=angle(ThetaHilbert);
ThetaAngle(~logicalValidLFPs)=NaN;

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAmp(~logicalValidLFPs)=NaN;
DeltaAngle=angle(DeltaHilbert);
DeltaAngle(~logicalValidLFPs)=NaN;

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;
glmTbl.WellSpikes=movsum(spike_train,round(0.4*re_fs),"omitmissing")';
glmTbl.ThetaAmp=ThetaAmp';
% glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.ThetaAngle=ThetaAngle';
glmTbl.SinThetaAngle=sin(glmTbl.ThetaAngle);
glmTbl.CosThetaAngle=cos(glmTbl.ThetaAngle);
glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.DeltaAngle=DeltaAngle';
glmTbl.SinDeltaAngle=sin(glmTbl.DeltaAngle);
glmTbl.CosDeltaAngle=cos(glmTbl.DeltaAngle);

% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'SinThetaAngle + CosThetaAngle + SinDeltaAngle + CosDeltaAngle+'...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*'...
%     'SinThetaAngle*CosThetaAngle*SinDeltaAngle*CosDeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + '...
%     'ThetaAmp*ThetaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
modelspec='WellSpikes ~ ThetaAmp*(ThetaAngle+SinThetaAngle+CosThetaAngle)';
% modelspec=['WellSpikes ~ ThetaAmp*ThetaAngle*SinThetaAngle*CosThetaAngle' ...
%     '*DeltaAmp*DeltaAngle*SinDeltaAngle*CosDeltaAngle'];

mdl_zeroed=fitglm(glmTbl,modelspec,'Distribution','poisson');
%% Frequency Distribution


%% Burst duration histograms

spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs.mat");
spike_dyn=spike_dyn.well_spike_dynamics_table;

burstLength=spike_dyn.BurstDuration{spike_dyn.fi==5 & spike_dyn.channel_name=="E10"};

figure
histogram(burstLength,logspace(log10(5),log10(200),25))
ax=gca;
ax.XScale="log";
xticks(round(logspace(log10(5),log10(200),25)))
xlabel("Burst Duration (ms)")
ylabel("nBursts")
%% GLM AxonLFP-Well Spikes NaN outside theta regions Spikes Per Burst

% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/10)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

validLFPIndex=[];
for nEndPts=1:size(LFPEndPts,1)
    validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
end
logicalValidLFPs=zeros(1,length(re_t));
logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
logicalValidLFPs=logical(logicalValidLFPs);

ThetaHilbert=hilbert(theta);
ThetaAmp=abs(ThetaHilbert);
ThetaAmp(~logicalValidLFPs)=NaN;
ThetaAngle=angle(ThetaHilbert);
ThetaAngle(~logicalValidLFPs)=NaN;

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
DeltaAmp(~logicalValidLFPs)=NaN;
DeltaAngle=angle(DeltaHilbert);
DeltaAngle(~logicalValidLFPs)=NaN;

angle_edges=[-pi:pi/180*18:pi];
spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\E10_spikes.mat");
spikes=spikes.index/1000;
spike_train=zeros(1,length(t));
spike_train(ismembertol(t,spikes))=1;
glmTbl=table();
spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
spike_train=zeros(1,300*re_fs);
spike_train(spike_idx)=1;
glmTbl.WellSpikes=movsum(spike_train,round(0.25*re_fs),"omitmissing")';
glmTbl.ThetaAmp=ThetaAmp';
% glmTbl.ThetaAngle=discretize(resample(ThetaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.ThetaAngle=ThetaAngle';
glmTbl.SinThetaAngle=sin(glmTbl.ThetaAngle);
glmTbl.CosThetaAngle=cos(glmTbl.ThetaAngle);
glmTbl.DeltaAmp=DeltaAmp';
% glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
glmTbl.DeltaAngle=DeltaAngle';
glmTbl.SinDeltaAngle=sin(glmTbl.DeltaAngle);
glmTbl.CosDeltaAngle=cos(glmTbl.DeltaAngle);

% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'SinThetaAngle + CosThetaAngle + SinDeltaAngle + CosDeltaAngle+'...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*'...
%     'SinThetaAngle*CosThetaAngle*SinDeltaAngle*CosDeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + '...
%     'ThetaAmp*ThetaAngle'];
% modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
%     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
modelspec='WellSpikes ~ ThetaAmp*(ThetaAngle+SinThetaAngle+CosThetaAngle)';
% modelspec=['WellSpikes ~ ThetaAmp*ThetaAngle*SinThetaAngle*CosThetaAngle' ...
%     '*DeltaAmp*DeltaAngle*SinDeltaAngle*CosDeltaAngle'];

mdl_zeroed=fitglm(glmTbl,modelspec,'Distribution','normal');