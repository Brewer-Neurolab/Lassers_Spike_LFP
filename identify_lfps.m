function identify_lfps(data,fs,wave_nsamples,nsamples_combine_thresh)
%data should be a 1d vector
%fs should be the sampling rate, input as a scalar
%the min length of a wave should be specified in samples
%nsamples_combine_thresh is the minimum number of samples two waves have to
%be from one another, otherwise combine them

%get hilbert transform of the data
myHilbert=hilbert(data);
dataAmplitude=abs(myHilbert);

%perform a gaussian convolution on the amplitude of hilbert data to smooth
%evelope of the signal
% gw=gausswin(wave_nsamples)
% dataAmplitude=conv(dataAmplitude,gw,'same');

%perform thresholding
ampSTD=std(dataAmplitude);
lowThresh=1.5*ampSTD;
highThresh=2.5*ampSTD;
threshedAmps=dataAmplitude>


end