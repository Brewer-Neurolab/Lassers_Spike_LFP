function sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn,nYbin,thresh_mult,parent_dir)

thetaAmpThresh=std(LFPAmplitude)*thresh_mult;
ampEdges=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);
ampCenters=convert_edges_2_centers(log10(ampEdges));
ampCenters=10.^ampCenters;
angleEdges=[0:18:360];
angleCenters=convert_edges_2_centers(angleEdges);
angleEdges2Cycle=[-360:18:360];
logicalValidLFPs=logical(logicalValidLFPs);

figure
tf=tiledlayout(1,3);
% tf=tiledlayout(2,3);
nexttile
histogram(LFPAngles(logicalValidLFPs),angleEdges)
title("High Angle Distribution")
xlim([min(angleEdges),max(angleEdges)])
xticks(angleEdges)

nexttile
histogram(LFPAmplitude(logicalValidLFPs),ampEdges)
title("High Amplitude Distribution")
set(gca,"XScale","log")
xlim([min(ampEdges),max(ampEdges)])
xticks(ampEdges(1:2:end))
xticklabels(round(ampEdges(1:2:end)))

nexttile
scatter(LFPAngles(logicalValidLFPs),LFPAmplitude(logicalValidLFPs))
title("High Amplitude vs Angle")
ylim([min(ampEdges),max(ampEdges)])
yticks(ampEdges(1:2:end))
yticklabels(round(ampEdges(1:2:end)))
xlim([min(angleEdges),max(angleEdges)])
xticks(angleEdges)
set(gca,"YScale","log")

title(tf,sourceElec+" axon distributions")

nIter=100;
for nElec=18%1:length(targetElecs)
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
    logicalValidSpikes=logical(logicalValidSpikes);

    nBurstsCounter(nElec)=sum(logicalBurstStarts & logicalValidLFPs);

    % figure
    wellBurstAngles=LFPAngles(logicalValidSpikes & logicalValidLFPs);
    % wellBurstAngles=[wellBurstAngles-360,wellBurstAngles];
    wellBurstAmp=LFPAmplitude(logicalValidSpikes & logicalValidLFPs);

    wellBurstStartAngles=LFPAngles(logicalBurstStarts & logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & logicalValidLFPs);

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_spb=well_spb(burstIdx);

    % for nBursts=1:length(wellBurstStartAngles)
    %     repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
    %     repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
    % end

    figure
    tf=tiledlayout(3,2,"Padding","tight","TileSpacing","tight");
    % well spikes vs amplitude
    nexttile
    histogram(wellBurstAmp,ampEdges)
    ampProbs=histcounts(wellBurstStartAmp,ampEdges,"Normalization","probability");
    ampMI=modulationIndex(ampProbs);
    ampPval=shuffleLFP_ModIdx(ampMI,LFPAmplitude,find(logicalValidLFPs),logicalValidSpikes,ampEdges,nIter);
    set(gca,"XScale","log")
    % axis square
    pbaspect([2,1,1])
    xlim([min(ampEdges),max(ampEdges)])
    xticks(ampEdges(1:2:end))
    xticklabels(round(ampEdges(1:2:end)))
    title("Well Spikes VS Amplitude")
    if ampPval<1/nIter
        subtitle("mod idx="+round(ampMI,2)+" p<"+1/nIter)
    else
        subtitle("mod idx="+round(ampMI,2)+" p="+ampPval)
    end
    ylabel("Spikes")
    xlabel("Amplitude uV")
    set(gca,"FontSize",12)

    nexttile
    binnedSPBAmp=discretize(wellBurstStartAmp,ampEdges);
    spbAmpMean=[];
    spbAmpSE=[];
    for nBins=1:length(ampCenters)
        if sum(binnedSPBAmp==nBins)>0
            spbAmpMean(nBins)=mean(well_spb(binnedSPBAmp==nBins));
            spbSD=std(well_spb(binnedSPBAmp==nBins));
            spbAmpSE(nBins)=spbSD/sqrt(sum(binnedSPBAmp==nBins));
        else
            spbAmpMean(nBins)=0;
            spbAmpSE(nBins)=0;
        end
    end
    % bar(ampCenters,spbAmpMean,1)
    h=histogram();
    h.BinEdges=ampEdges;
    h.BinCounts=spbAmpMean;
    hold on
    errorbar(ampCenters,spbAmpMean,spbAmpSE,"LineStyle","none","Color",'k')
    % scatter(wellBurstStartAmp,well_spb(burstIdx))
    set(gca,"XScale","log")
    % axis square
    pbaspect([2,1,1])
    xlim([min(ampEdges),max(ampEdges)])
    xticks(ampEdges(1:2:end))
    xticklabels(round(ampEdges(1:2:end)))
    title("Well Spikes Per Burst VS Burst Start LFP Amplitude")
    
    xlabel("Amplitude uV")
    set(gca,"FontSize",12)
    ylabel("Spikes per Burst","FontSize",10)
    currentYLim=ylim;
    ylim([min(currentYLim),max(currentYLim)*1.1])

    %well spikes vs angle
    nexttile
    histogram([wellBurstAngles-360,wellBurstAngles],angleEdges2Cycle)
    angleProbs=histcounts(wellBurstStartAngles,angleEdges,"Normalization","probability");
    angleMI=modulationIndex(angleProbs);
    anglePval=shuffleLFP_ModIdx(angleMI,LFPAngles,find(logicalValidLFPs),logicalValidSpikes,angleEdges,nIter);
    % axis square
    pbaspect([2,1,1])
    xlim([-180,180])
    xticks(angleEdges2Cycle(1:2:end))
    title("Well Spikes VS Angle")
    if anglePval<1/nIter
        subtitle("mod idx="+round(angleMI,2)+" p<"+1/nIter)
    else
        subtitle("mod idx="+round(angleMI,2)+" p="+ampPval)
    end
    ylabel("Spikes")
    xlabel("Angle")
    set(gca,"FontSize",12)

    nexttile
    binnedSPBAngle=discretize(wellBurstStartAngles,angleEdges);
    spbAngleMean=[];
    spbAngleSE=[];
    for nBins=1:length(ampCenters)
        if sum(binnedSPBAngle==nBins)>0
            spbAngleMean(nBins)=mean(well_spb(binnedSPBAngle==nBins));
            spbSD=std(well_spb(binnedSPBAngle==nBins));
            spbAngleSE(nBins)=spbSD/sqrt(sum(binnedSPBAngle==nBins));
        else
            spbAngleMean(nBins)=0;
            spbAngleSE(nBins)=0;
        end
    end
    h=histogram();
    h.BinEdges=angleEdges;
    h.BinCounts=spbAngleMean;
    hold on
    errorbar(angleCenters,spbAngleMean,spbAngleSE,"LineStyle","none","Color",'k')
    % scatter(wellBurstStartAngles,well_spb(burstIdx))
    % axis square
    pbaspect([2,1,1])
    xlim([min(angleEdges),max(angleEdges)])
    xticks(angleEdges)
    title("Well Spikes Per Burst VS Burst Start LFP Angle")
    
    xlabel("Angles")
    set(gca,"FontSize",12)
    ylabel("Spikes per Burst","FontSize",10)

    %axon amplitude and angle at well spikes
    % nexttile([1,2])
    % scatter(wellBurstAmp,wellBurstAngles)
    % set(gca,"XScale","log")
    % axis square
    % xlim([min(ampEdges),max(ampEdges)])
    % xticks(ampEdges(1:2:end))
    % xticklabels(round(ampEdges(1:2:end)))
    % ylim([min(angleEdges),max(angleEdges)])
    % yticks(angleEdges(1:2:end))
    % title("Well Spike Angle vs Well Amplitude")
    % ylabel("Angle")
    % xlabel("Amplitude uV")
    % set(gca,"FontSize",16)
    X=[[wellBurstAngles,-wellBurstAngles];[wellBurstAmp,wellBurstAmp]]';
    edges={[-360:18:360],logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)};
    if ~isempty(X)
        %calculates pxy
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells_xy{nElec}=N(1:40,1:nYbin);
        binxcenters{nElec}=convert_edges_2_centers([-360:18:360]);
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)));
        binxedges{nElec}=[-360:18:360];
        binyedges{nElec}=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);

        %calculate px
        bincount_cells_x{nElec}=histcounts(wellBurstAngles,[-360:18:360]);

        %calculate py
        bincount_cells_y{nElec}=histcounts(wellBurstAmp,logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1));

        % figure
        % histogram(wellBurstAngles,[0:18:360],"Normalization","probability")
        % xlabel("Angle degrees")
        % ylabel("Weighted Burst Probability in Soma")
        % title(targetElecs(nElec)+" n="+length(wellBurstAmp))
        % axis square
        % xticks([0:36:360])
        % xlim([0,360])
        % set(gca,"FontSize",16)
        % 
        % figure
        % histogram(wellBurstAmp,logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1),"Normalization","probability")
        % xlabel("Amplitude uV")
        % ylabel("Weighted Burst Probability in Soma")
        % ax=gca;
        % ax.XScale="log";
        % title(targetElecs(nElec)+" n="+length(wellBurstAmp))
        % axis square
        % xticks(round(logspace(1,4,31)))
        % % xlim(1,10000)
        % set(gca,"FontSize",16)

    else
        bincount_cells_xy{nElec}=[];
        binxcenters{nElec}=[];
        binycenters{nElec}=[];
        binxedges{nElec}=[];
        binyedges{nElec}=[];
    end

    % nexttile([3,4])
    nexttile
    maxHeightAll=cellfun(@(x) max(x,[],"all"),bincount_cells_xy,'UniformOutput',false);
    maxHeightAll(cellfun(@isempty,maxHeightAll))={0};
    maxHeightAll=max(cell2mat(maxHeightAll));

    myCLim=[0,maxHeightAll];
    imagesc(binxcenters{nElec},binycenters{nElec},flipud(rot90(bincount_cells_xy{nElec},1)),myCLim)
    ylim([min(binyedges{nElec}),max(binyedges{nElec})])
    xlim([-180,180])
    % axis square
    pbaspect([2,1,1])

    set(gca,"FontSize",18)
    xticks(binxedges{nElec}(1:2:end))
    yticks(binyedges{nElec}(1:4:end))
    yticklabels(round(binyedges{nElec}(1:4:end),2,"significant"))
    set(gca,'YDir','normal')
    set(gca,"YScale","log")

    cb = colorbar("FontSize",18);
    colormap hot
    % cb.Layout.Tile = 'east';
    cb.Limits=myCLim;
    clim(myCLim)
    cb.Ticks=round(linspace(0,maxHeightAll,5));

    title(tf,targetElecs(nElec),"FontSize",18)
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
end

end