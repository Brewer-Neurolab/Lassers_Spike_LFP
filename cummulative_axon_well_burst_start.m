function cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn)

bincount_cells=[];
binxcenters=[];
binxedges=[];
binycenters=[];
binyedges=[];

for nElec=1:length(targetElecs)
    disp("Calculating "+targetElecs(nElec))
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_burst_starts=well_burst_bounds(:,1);
    % remap burst starts to new sampling rate
    well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
    logicalBurstStarts=zeros(1,length(re_t));
    logicalBurstStarts(well_burst_starts)=1;

    % figure
    wellBurstStartAngles=LFPAngles(logicalBurstStarts & logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & logicalValidLFPs);

    thetaAmpThresh=std(LFPAmplitude);

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};

    for nBursts=1:length(wellBurstStartAngles)
        repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
        repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
    end

    X=[repwellBurstStartAngles;repwellBurstStartAmp]';
    edges={[0:40:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10)};

    if ~isempty(X)
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells{nElec}=N(1:9,1:9);
        binxcenters{nElec}=convert_edges_2_centers([0:40:360]);
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10)));
        binxedges{nElec}=[0:40:360];
        binyedges{nElec}=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10);
    else
        bincount_cells{nElec}=[];
        binxcenters{nElec}=[];
        binycenters{nElec}=[];
        binxedges{nElec}=[];
        binyedges{nElec}=[];
    end

end

phase_amp_heatmap_formatter(bincount_cells,binxcenters,binxedges,binycenters,binyedges,sourceElec,targetElecs)

end