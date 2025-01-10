function cummulative_axon_well_burst_start(LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn)

bincount_cells=[];
binx=[];
biny=[];

for nElec=1:length(targetElecs)
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_burst_starts=well_burst_bounds(:,1);
    % remap burst starts to new sampling rate
    well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
    logicalBurstStarts=zeros(1,length(re_t));
    logicalBurstStarts(well_burst_starts)=1;

    figure
    wellBurstStartAngles=LFPAngles(logicalBurstStarts & logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & logicalValidLFPs);

    thetaAmpThresh=std(LFPAmplitude);

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==6 & well_spike_dyn.channel_name=="E10"};

    for nBursts=1:length(wellBurstStartAngles)
        repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
        repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
    end

    X=[repwellBurstStartAngles;repwellBurstStartAmp]';
    edges={[0:40:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),10)};

    [N,C]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
    % xlabel("Axon Phase Angle (deg)")
    % ylabel("Axon Burst Start Amplitude (uV)")
    % zlabel("Cummulative Soma Spikes Per Burst")
    % 
    % xlim([0,360])
    % yticks(logspace(log10(min(repwellBurstStartAmp,[],"all")),log10(max(repwellBurstStartAmp,[],"all")),10))
    % yticklabels(round(logspace(log10(min(repwellBurstStartAmp,[],"all")),log10(max(repwellBurstStartAmp,[],"all")),10),1))
    % ylim([min(repwellBurstStartAmp,[],"all"),max(repwellBurstStartAmp,[],"all")])
    % 
    % ax=gca;
    % ax.YScale="log";
    bincount_cells{nElec}=N;
    binx{nElec}=convert_edges_2_centers([0:40:360]);
    biny{nElec}=10.^convert_edges_2_centers(logspace(log10(min(repwellBurstStartAmp,[],"all")),log10(max(repwellBurstStartAmp,[],"all")),10));
end

phase_amp_heatmap_formatter(bincount_cells,binx,biny,sourceElec,targetElecs)

end