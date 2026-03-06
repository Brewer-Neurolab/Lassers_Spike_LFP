function [pval,MIVec,myProbability]=fftShuffleLFP_ModIdx(myMI,myData,LFPEndPts,spikeIdx,myEdges,nIter,hilbertDim)

% define regions to shuffle from high amplitude true idxs

dataCells=[];
spikeCells=[];
for nEndPts=1:size(LFPEndPts,1)
    dataCells{nEndPts}=myData(LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2));
    spikeCells{nEndPts}=spikeIdx(LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2));
    % spikeCells=[spikeCells,spikeIdx(LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2))];
end

% spikeCells=logical(spikeCells);

MIVec=[];
myProbability=[];
% https://www.mathworks.com/matlabcentral/answers/451578-fft-and-ifft-random-phases
for myIter=1:nIter
    randData=[];
    % rand_theta=2*rand(1,1);
    for nSegments=1:length(dataCells)
        % fft angle shift
        % n=length(dataCells{nSegments});

        % f=fft(dataCells{nSegments});
        % complexf=1i*f;
        % complexf=f;

        % posF=2:floor(n/2)+mod(n,2);
        % negF=ceil(n/2)+1+~mod(n,2):n;

        % if mod(length(dataCells{nSegments}),2)==1
        %     rand_theta=2*rand(1,floor(length(dataCells{nSegments})/2));
        % else
        %     rand_theta=2*rand(1,floor(length(dataCells{nSegments})/2-1));
        % end

        % Multiplying phases by random integer is insufficient
        % f(posF)=f(posF) + -rand_theta.*1i.*complexf(posF);
        % f(negF)=f(negF) + fliplr(rand_theta).*1i.*complexf(negF);
        % f(negF)=f(negF) + rand_theta.*1i.*complexf(negF);

        % get imaginary component and randomize

        % f_rotate=ifft(f,"symmetric");

        f_rotate=phaseran(dataCells{nSegments}',1)';
        spikeCells{nSegments}=spikeCells{nSegments}(1:length(f_rotate));

        if hilbertDim=="amp"
            randData=[randData,abs(hilbert(f_rotate))];
        elseif hilbertDim=="angle"
            randData=[randData,wrapTo180(rad2deg(angle(hilbert(f_rotate))))];
        else
            randData=[randData,f_rotate];
        end

    end

    spikelocs=cell2mat(spikeCells);

    myProbability(nIter,:)=histcounts(randData(spikelocs),myEdges,"Normalization","probability");
    MI=modulationIndex(myProbability(nIter,:));

    MIVec(myIter)=MI;
end

if ~isempty(MIVec)
    pval=sum(MIVec>=myMI)/length(MIVec);
else
    pval=NaN;
end

end