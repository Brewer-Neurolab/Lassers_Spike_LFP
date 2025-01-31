function miTable=cummulative_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn,nYbin)

bincount_cells_xy=[];
bincount_cells_x=[];
bincount_cells_y=[];
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

    nBurstsCounter(nElec)=sum(logicalBurstStarts & logicalValidLFPs);

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
    edges={[0:40:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)};

    if ~isempty(X)
        %calculates pxy
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells_xy{nElec}=N(1:nYbin,1:9);
        binxcenters{nElec}=convert_edges_2_centers([0:40:360]);
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)));
        binxedges{nElec}=[0:40:360];
        binyedges{nElec}=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);

        %calculate px
        bincount_cells_x{nElec}=histcounts(repwellBurstStartAngles,[0:40:360]);

        %calculate py
        bincount_cells_y{nElec}=histcounts(repwellBurstStartAmp,logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1));
    else
        bincount_cells_xy{nElec}=[];
        binxcenters{nElec}=[];
        binycenters{nElec}=[];
        binxedges{nElec}=[];
        binyedges{nElec}=[];
    end

end

sourceProps.LFPAmps=LFPAmplitude;
sourceProps.LFPAngles=LFPAngles;
sourceProps.logicalValidLFPs=logicalValidLFPs;
sourceProps.well_spike_dyn=well_spike_dyn;
sourceProps.fi=fi;
sourceProps.t=t;
sourceProps.re_t=re_t;
sourceProps.nIter=100;

%comment any out for testing
% miTable=phase_amp_heatmap_formatter(bincount_cells_xy,binxcenters,binxedges,binycenters,binyedges,sourceElec,targetElecs,sourceProps);
% miTable=phase_amp_heatmap_formatter_percent(bincount_cells,binxcenters,binxedges,binycenters,binyedges,sourceElec,targetElecs,sourceProps);
miTable=mutualInfo_heatmap_formatter(bincount_cells_xy,binxcenters,binxedges,...
    binycenters,binyedges,bincount_cells_x,bincount_cells_y,sourceElec,targetElecs,sourceProps,nBurstsCounter);

end