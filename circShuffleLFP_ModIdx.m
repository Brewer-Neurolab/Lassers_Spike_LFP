function [pval,MIVec]=circShuffleLFP_ModIdx(myMI,myData,LFPEndPts,spikeIdx,myEdges,nIter)

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
% https://www.mathworks.com/matlabcentral/answers/451578-fft-and-ifft-random-phases
for myIter=1:nIter
    randData=[];
    rand_shift=randi(length(myData));
    for nSegments=1:length(dataCells)

        dat_rotate=circshift(dataCells{nSegments}',rand_shift)';
        spikeCells{nSegments}=spikeCells{nSegments}(1:length(dat_rotate));

        randData=[randData,dat_rotate];

    end

    spikelocs=cell2mat(spikeCells);

    myProbability=histcounts(randData(spikelocs),myEdges,"Normalization","probability");
    MI=modulationIndex(myProbability);

    MIVec(myIter)=MI;
end

if ~isempty(MIVec)
    pval=sum(MIVec>=myMI)/length(MIVec);
else
    pval=NaN;
end

end