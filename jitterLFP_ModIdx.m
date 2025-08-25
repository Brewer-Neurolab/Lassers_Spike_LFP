function [pval,MIVec]=jitterLFP_ModIdx(myMI,myData,dataIdx,spikeIdx,myEdges,nIter,jitterMs,fs)

%shuffles the modulation index
rng('default')
MIVec=[];

%convert jitter ms to samples
jitter=round(jitterMs/1000*fs);

for i=1:nIter

    %jitter by randomizing spike idx
    jitterVec=randi([-jitter,jitter],[sum(spikeIdx),1]);
    randData=myData(ismembertol(dataIdx,find(spikeIdx')+jitterVec,10^-12));
    
    myProbability=histcounts(randData,myEdges,"Normalization","probability");
    
    %calculate modulation index
    MI=modulationIndex(myProbability);

    MIVec(i)=MI;

end

if ~isempty(MIVec)
    pval=sum(MIVec>=myMI)/length(MIVec);
else
    pval=NaN;
end

end