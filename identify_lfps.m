function [endPts,dataAmplitude,myHilbert]=identify_lfps(data,fs, t_s,thresh,minLFPLength,minLFPCycles,nsamples_combine_thresh)
%data should be a 1d vector
%fs should be the sampling rate, input as a scalar
%the min length of a wave should be specified in samples
%nsamples_combine_thresh is the minimum number of samples two waves have to
%be from one another, otherwise combine them

%define time vector for plotting
t=0:1/fs:t_s-(1/fs);

%get hilbert transform of the data
myHilbert=hilbert(data);
dataAmplitude=abs(myHilbert);

%perform a gaussian convolution on the amplitude of hilbert data to smooth
%evelope of the signal
gw=gausswin(round(minLFPLength/minLFPCycles),5);
dataAmplitude=conv(dataAmplitude,gw,'same');

%rescale data amplitude
dataAmplitude=dataAmplitude*(max(abs(myHilbert))/max(dataAmplitude));

%perform thresholding
ampSTD=std(dataAmplitude);
lowThresh=thresh*ampSTD;
% highThresh=2.5*ampSTD;
isAboveThresh=dataAmplitude>=lowThresh;

% combine LFPs
nChanges=Inf;
while nChanges~=0
    nChanges=0;
    startstopidx=find(diff(isAboveThresh)~=0);

    if mod(length(startstopidx),2)~=0
        if isAboveThresh(1)
            startstopidx(1)=[];
        elseif isAboveThresh(end)
            startstopidx(end)=[];
        end
    end

    endPts=reshape(startstopidx,2,[])';
    for nGaps=1:size(endPts,1)-1
        gapSize=endPts(nGaps+1,1)-endPts(nGaps,2);
        if gapSize<=nsamples_combine_thresh
            isAboveThresh(endPts(nGaps,2):endPts(nGaps+1,1))=1;
            nChanges=nChanges+1;
        end
    end
end

%eliminate LFPs that are too short
lfpLen=diff(endPts,1,2);
endPts(lfpLen<minLFPLength,:)=[];

figure
hold on
plot(t,data)
plot(t,dataAmplitude,'k')
plot(t(endPts(:,1)),dataAmplitude(endPts(:,1)),'g.',"MarkerSize",20)
plot(t(endPts(:,2)),dataAmplitude(endPts(:,2)),'r.',"MarkerSize",20)
hold off
end