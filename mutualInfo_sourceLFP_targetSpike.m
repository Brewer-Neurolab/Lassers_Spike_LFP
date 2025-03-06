function miTable=mutualInfo_sourceLFP_targetSpike(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn,nYbin,thresh_mult,parent_dir)

bincount_cells_xy=[];
bincount_cells_x=[];
bincount_cells_y=[];
binxcenters=[];
binxedges=[];
binycenters=[];
binyedges=[];

for nElec=1:length(targetElecs)
    well_spikes=load(parent_dir+targetElecs(nElec)+"_spikes.mat");
    well_spikes=well_spikes.index/1000; % in seconds
    well_spikes_idx=find(ismembertol(t,well_spikes));
    % well_spikes=remap(well_spikes,1/1000,300,1/25000,300);
    
    disp("Calculating "+targetElecs(nElec))
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_burst_starts=well_burst_bounds(:,1);
    % remap burst starts to new sampling rate
    well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
    logicalBurstStarts=zeros(1,length(re_t));
    logicalBurstStarts(well_burst_starts)=1;
    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    myBurstBounds=well_burst_bounds(burstIdx,:);
    validWellSpikes=[];
    for nBursts=1:size(myBurstBounds,1)
        validWellSpikes=[validWellSpikes,well_spikes_idx(well_spikes_idx>=myBurstBounds(nBursts,1) & well_spikes_idx<=myBurstBounds(nBursts,2))];
    end
    % well_spikes_idx=t(validWellSpikes);
    validWellSpikes=round(remap(validWellSpikes,1,length(t),1,length(re_t)));
    logicalValidSpikes=zeros(1,length(re_t));
    logicalValidSpikes(validWellSpikes)=1;

    nBurstsCounter(nElec)=sum(logicalBurstStarts & logicalValidLFPs);

    % figure
    wellBurstAngles=LFPAngles(logicalValidSpikes & logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstAmp=LFPAmplitude(logicalValidSpikes & logicalValidLFPs);

    thetaAmpThresh=std(LFPAmplitude)*thresh_mult;

    % %repeat for spikes per burst
    % repwellBurstStartAngles=[];
    % repwellBurstStartAmp=[];
    % 
    
    % well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    % 
    % for nBursts=1:length(wellBurstAngles)
    %     repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
    %     repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
    % end

    X=[wellBurstAngles;wellBurstAmp]';
    edges={[0:18:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)};

    if ~isempty(X)
        %calculates pxy
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells_xy{nElec}=N(1:20,1:nYbin);
        binxcenters{nElec}=convert_edges_2_centers([0:18:360]);
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)));
        binxedges{nElec}=[0:18:360];
        binyedges{nElec}=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);

        %calculate px
        bincount_cells_x{nElec}=histcounts(wellBurstAngles,[0:18:360]);

        %calculate py
        bincount_cells_y{nElec}=histcounts(wellBurstAmp,logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1));

        figure
        histogram(wellBurstAngles,[0:18:360],"Normalization","probability")
        xlabel("Angle degrees")
        ylabel("Weighted Burst Probability in Soma")
        title(targetElecs(nElec)+" n="+length(wellBurstAmp))
        axis square
        xticks([0:36:360])
        xlim([0,360])
        set(gca,"FontSize",16)

        figure
        histogram(wellBurstAmp,logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1),"Normalization","probability")
        xlabel("Amplitude uV")
        ylabel("Weighted Burst Probability in Soma")
        ax=gca;
        ax.XScale="log";
        title(targetElecs(nElec)+" n="+length(wellBurstAmp))
        axis square
        xticks(round(logspace(1,4,31)))
        % xlim(1,10000)
        set(gca,"FontSize",16)

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
    binycenters,binyedges,bincount_cells_x,bincount_cells_y,sourceElec,targetElecs,sourceProps,nBurstsCounter,nYbin);

end