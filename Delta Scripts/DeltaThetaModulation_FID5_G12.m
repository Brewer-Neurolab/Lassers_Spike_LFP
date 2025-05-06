clear
clc
close all
%% Axon Setup
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
theta=filtfilt(sos,g,LFP_filt);

%filter for delta
[A,B,C,D]=butter(8,4/(re_fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,0.5/(re_fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
delta=filtfilt(sos,g,LFP_filt);

%get axon spikes
axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
axon_spikes=axon_spikes.cluster_class(:,2)/1000;
axon_spike_train=zeros(1,length(t));
axon_spike_train(ismembertol(t,axon_spikes))=1;

nperm=100;
%% Get params for modulation Index

% plot axon data tagged
thresh_mult=1;

%define max number of samples for combining LFPs
nsamples_combine_thresh=(1/4)*re_fs*3; 
% nsamples_combine_thresh=[];

%define min lfp length as 2x shortest theta cycle
minLFPCycles=2; %default 2
minLFPLength=(1/10)*minLFPCycles*re_fs;

[LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(delta,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
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
% ThetaAmp(~logicalValidLFPs)=NaN;
ThetaAngle=angle(ThetaHilbert);
% ThetaAngle(~logicalValidLFPs)=NaN;

DeltaHilbert=hilbert(delta);
DeltaAmp=abs(DeltaHilbert);
% DeltaAmp(~logicalValidLFPs)=NaN;
DeltaAngle=angle(DeltaHilbert);
% DeltaAngle(~logicalValidLFPs)=NaN;

%Amp-Amp MI
nThetaAmpBins=freedmandiaconis(ThetaAmp);
nDeltaAmpBins=freedmandiaconis(DeltaAmp);
nThetaAngleBins=freedmandiaconis(ThetaAngle);
nDeltaAngleBins=freedmandiaconis(DeltaAngle);

% 40 optimal bins

thetaAmpEdges=logspace(log10(min(ThetaAmp)),log10(max(ThetaAmp)),41);
deltaAmpEdges=logspace(log10(min(DeltaAmp)),log10(max(DeltaAmp)),41);
thetaAngleEdges=linspace(-pi,pi,41);
deltaAngleEdges=linspace(-pi,pi,41);

%% calculate delta phase theta amp MI

discretizeDelta=discretize(DeltaAngle,deltaAngleEdges);

pDAngleTAmp=zeros([1,40]);
for i=1:40
    pDAngleTAmp(i)=mean(ThetaAmp(discretizeDelta==i),"omitmissing");
end
pDAngleTAmp=pDAngleTAmp/sum(pDAngleTAmp);

figure
histogram('BinEdges',linspace(-180,180,41),'BinCounts',pDAngleTAmp)
xlabel("Delta Angle")
ylabel("Theta Amplitude")
xticks(-180:18:180)
set(gca,"FontSize",20)

MIDAngleTAmp=modulationIndex(pDAngleTAmp);

% Shuffle theta amp to determine significance of MI
miVec=[];
for n=1:nperm
    shuffleThetaAmp=ThetaAmp(randperm(length(ThetaAmp)));
    pDAngleTAmp=zeros([1,40]);
    for i=1:40
        pDAngleTAmp(i)=mean(shuffleThetaAmp(discretizeDelta==i),"omitmissing");
    end
    pDAngleTAmp=pDAngleTAmp/sum(pDAngleTAmp);
    miVec(n)=modulationIndex(pDAngleTAmp);
end

pMIDAngleTAmp=sum(miVec>=MIDAngleTAmp)/nperm;

%% calculate delta amp theta amp MI

discretizeDelta=discretize(DeltaAmp,deltaAmpEdges);

pDAmpTAmp=zeros([1,40]);
for i=1:40
    pDAmpTAmp(i)=mean(ThetaAmp(discretizeDelta==i),"omitmissing");
end
pDAmpTAmp=pDAmpTAmp/sum(pDAmpTAmp);

figure
histogram('BinEdges',deltaAmpEdges,'BinCounts',pDAmpTAmp)
xlabel("Delta Amplitude")
ylabel("Theta Amplitude")
set(gca,"FontSize",20)
xscale(gca,"log")
xlim([min(deltaAmpEdges),max(deltaAmpEdges)])

MIDAmpTAmp=modulationIndex(pDAmpTAmp);

% Shuffle theta amp to determine significance of MI
miVec=[];
for n=1:nperm
    shuffleThetaAmp=ThetaAmp(randperm(length(ThetaAmp)));
    pDAmpTAmp=zeros([1,40]);
    for i=1:40
        pDAmpTAmp(i)=mean(shuffleThetaAmp(discretizeDelta==i),"omitmissing");
    end
    pDAmpTAmp=pDAmpTAmp/sum(pDAmpTAmp);
    miVec(n)=modulationIndex(pDAmpTAmp);
end

pMIDAmpTAmp=sum(miVec>=MIDAmpTAmp)/nperm;
%% calculate theta phase delta amp MI

discretizeTheta=discretize(ThetaAngle,thetaAngleEdges);

pTAngleDAmp=zeros([1,40]);
for i=1:40
    pTAngleDAmp(i)=mean(DeltaAmp(discretizeTheta==i),"omitmissing");
end
pTAngleDAmp=pTAngleDAmp/sum(pTAngleDAmp);

figure
histogram('BinEdges',linspace(-180,180,41),'BinCounts',pTAngleDAmp)
xlabel("Theta Angle")
ylabel("Delta Amplitude")
xticks(-180:18:180)
set(gca,"FontSize",20)

MITAngleDAmp=modulationIndex(pTAngleDAmp);

% Shuffle theta amp to determine significance of MI
miVec=[];
for n=1:nperm
    shuffleDeltaAmp=DeltaAmp(randperm(length(DeltaAmp)));
    pTAngleDAmp=zeros([1,40]);
    for i=1:40
        pTAngleDAmp(i)=mean(shuffleDeltaAmp(discretizeTheta==i),"omitmissing");
    end
    pTAngleDAmp=pTAngleDAmp/sum(pTAngleDAmp);
    miVec(n)=modulationIndex(pTAngleDAmp);
end

pMITAngleDAmp=sum(miVec>=MIDAngleTAmp)/nperm;

%% calculate delta amp theta amp MI

discretizeTheta=discretize(ThetaAmp,thetaAmpEdges);

pTAmpDAmp=zeros([1,40]);
for i=1:40
    if ~isnan(mean(DeltaAmp(discretizeTheta==i),"omitmissing"))
        pTAmpDAmp(i)=mean(DeltaAmp(discretizeTheta==i),"omitmissing");
    end
end
pTAmpDAmp=pTAmpDAmp/sum(pTAmpDAmp);

figure
histogram('BinEdges',thetaAmpEdges,'BinCounts',pTAmpDAmp)
xlabel("Theta Amplitude")
ylabel("Delta Amplitude")
set(gca,"FontSize",20)
xscale(gca,"log")
xlim([min(thetaAmpEdges),max(thetaAmpEdges)])

MITAmpDAmp=modulationIndex(pTAmpDAmp);

% Shuffle theta amp to determine significance of MI
miVec=[];
for n=1:nperm
    shuffleDeltaAmp=DeltaAmp(randperm(length(DeltaAmp)));
    pTAmpDAmp=zeros([1,40]);
    for i=1:40
        pTAmpDAmp(i)=mean(shuffleDeltaAmp(discretizeTheta==i),"omitmissing");
    end
    pTAmpDAmp=pTAmpDAmp/sum(pTAmpDAmp);
    miVec(n)=modulationIndex(pTAmpDAmp);
end

pMITAmpDAmp=sum(miVec>=MITAmpDAmp)/nperm;