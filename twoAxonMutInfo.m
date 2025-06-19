function [MutInfo,pval]=twoAxonMutInfo(dir1,dir2,HPHz,LPHz,fs,re_fs,t_rec,nperm)
pval=[];
%% Axon 1 Setup
data=load(dir1);
data=data.data;

re_t=0:1/re_fs:t_rec-(1/re_fs);

%downsample data
data=resample(data,re_fs,fs);

%create full sampled time steps
t=0:1/fs:t_rec-(1/fs);

%filter for theta
[A,B,C,D]=butter(8,LPHz/(re_fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
LFP_filt=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,HPHz/(re_fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
axon1_LFP=filtfilt(sos,g,LFP_filt);

%% Axon 2 Axon Setup
data=load(dir2);
data=data.data;

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
axon2_LFP=filtfilt(sos,g,LFP_filt);

%% Mutual info Axon1 Axon2 FID 5

Axon1Hilbert=hilbert(axon1_LFP);
Axon1Amp=abs(Axon1Hilbert);
Axon1Angle=angle(Axon1Hilbert);
Axon2Hilbert=hilbert(axon2_LFP);
Axon2Amp=abs(Axon2Hilbert);
Axon2Angle=angle(Axon2Hilbert);

% Get optimal number of bins TODO

n_bins=20;

figure
edgesAxon1=logspace(log10(min(Axon1Amp)),log10(max(Axon1Amp)),n_bins+1);
hAxon1=histogram(Axon1Amp,edgesAxon1);
probAxon1=hAxon1.Values./sum(hAxon1.Values);
SEAxon1=-sum(probAxon1.*log2(probAxon1+eps));
xscale(gca,"log")
xlim([min(Axon1Amp),max(Axon1Amp)])

figure
edgesAxon2=logspace(log10(min(Axon2Amp)),log10(max(Axon2Amp)),n_bins+1);
hAxon2=histogram(Axon2Amp,edgesAxon2);
probAxon2=hAxon2.Values./sum(hAxon2.Values);
SEAxon2=-sum(probAxon2.*log2(probAxon2+eps));
xscale(gca,"log")
xlim([min(Axon2Amp),max(Axon2Amp)])

[~,~,bins1]=histcounts(Axon1Amp,edgesAxon1);
[~,~,bins2]=histcounts(Axon2Amp,edgesAxon2);

jp=zeros(n_bins);
for i1=1:n_bins
    for i2=1:n_bins
        jp(i1,i2)=sum(bins1==i1 & bins2==i2);
    end
end

jp=jp./sum(sum(jp));

JSE=-sum(sum(jp.*log2(jp+eps)));

MutInfo=SEAxon1+SEAxon2-JSE;

%Nonparametric stats
figure
MIVec=[];
for i=1:nperm

    randAmp1=Axon1Amp(randperm(length(Axon1Amp)));
    edgesAxon1=logspace(log10(min(randAmp1)),log10(max(randAmp1)),n_bins+1);
    hAxon1=histogram(randAmp1,edgesAxon1);
    probAxon1=hAxon1.Values./sum(hAxon1.Values);
    SEAxon1=-sum(probAxon1.*log2(probAxon1+eps));
    xscale(gca,"log")
    xlim([min(randAmp1),max(randAmp1)])

    randAmp2=Axon2Amp;%(randperm(length(Axon2Amp)));
    edgesAxon2=logspace(log10(min(randAmp2)),log10(max(randAmp2)),n_bins+1);
    hAxon2=histogram(randAmp2,edgesAxon2);
    probAxon2=hAxon2.Values./sum(hAxon2.Values);
    SEAxon2=-sum(probAxon2.*log2(probAxon2+eps));
    xscale(gca,"log")
    xlim([min(Axon2Amp),max(Axon2Amp)])

    [~,~,bins1]=histcounts(randAmp1,edgesAxon1);
    [~,~,bins2]=histcounts(randAmp2,edgesAxon2);

    jp=zeros(n_bins);
    for i1=1:n_bins
        for i2=1:n_bins
            jp(i1,i2)=sum(bins1==i1 & bins2==i2);
        end
    end

    jp=jp./sum(sum(jp));

    JSE=-sum(sum(jp.*log2(jp+eps)));

    MIVec(i)=SEAxon1+SEAxon2-JSE;
end
close gcf

pval=sum(MIVec>=MutInfo)/nperm;

if pval>=0.05
    MutInfo=NaN;
end


end