function [pval,MI_Vec]=shuffleCMI(myMI,LFPAmps,LFPAngles,logicalValidLFPs,well_spike_dyn,myElec,fi,t,re_t,nIter,nYbin)
MI_Vec=[];
pval=[];

well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==myElec};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;
thetaAmpThresh=std(LFPAmps);

rng('default')
nPermutations=nIter;
shuffledLFPAngles=zeros(nPermutations,length(LFPAngles));
shuffledLFPAmps=zeros(nPermutations,length(LFPAmps));

%calculate MI for each permuatation
for nPerm=1:nPermutations
    randPermAngles=randperm(length(LFPAngles));
    shuffledLFPAngles(nPerm,:)=LFPAngles(randPermAngles);
    randPermAmps=randperm(length(LFPAmps));
    shuffledLFPAmps(nPerm,:)=LFPAmps(randPermAmps);

    wellBurstStartAngles=shuffledLFPAngles(nPerm,logicalBurstStarts & logicalValidLFPs);
    wellBurstStartAmp=shuffledLFPAmps(nPerm,logicalBurstStarts & logicalValidLFPs);

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==myElec};

    for nBursts=1:length(wellBurstStartAngles)
        repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
        repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
    end

    X=[repwellBurstStartAngles;repwellBurstStartAmp]';
    edges={[0:18:360],logspace(log10(min(repwellBurstStartAmp)),log10(max(LFPAmps)),nYbin)};

    if ~isempty(X)
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincounts_cells=N(1:20,1:nYbin);
        binxcenters=convert_edges_2_centers([0:18:360]);
        binycenters=10.^convert_edges_2_centers(log10(logspace(log10(min(repwellBurstStartAmp)),log10(max(LFPAmps)),nYbin+1)));
        binxedges=[0:18:360];
        binyedges=logspace(log10(min(repwellBurstStartAmp)),log10(max(LFPAmps)),nYbin+1);
        bincounts_x=histcounts(repwellBurstStartAngles);
        bincounts_y=histcounts(repwellBurstStartAmp);
    else
        bincounts_cells=[];
        binxcenters=[];
        binycenters=[];
        binxedges=[];
        binyedges=[];
    end

    MI_Vec(nPerm)=mutualInfo(bincounts_cells/sum(bincounts_cells,"all"),...
            bincounts_x/sum(bincounts_x),...
            bincounts_y/sum(bincounts_y),...
            sum(bincounts_cells,"all"));

    %comment out when not testing
    % figure
    % imagesc(binxcenters,binycenters,rot90(flipud(bincounts_cells),-1))
    % xticks(binxcenters)
    % yticks(binycenters)
    % yticklabels(round(binycenters,2,"significant"))
    % set(gca,'YDir','normal')
    % set(gca,"YScale","log")
    % axis("square")
    % cb = colorbar;
    % cb.Limits=[0,200];
    % colormap hot
    % clim([0,200])
    % 
end

%calculate total pval
if ~isempty(MI_Vec)
    pval=sum(MI_Vec>=myMI)/length(MI_Vec);
end

% for testing and validating pval
% figure
% histogram(MI_Vec)
% hold on
% xline(myMI)
% hold off

end