function relationTable=sourceLFP_targetSpike_BurstRelations(t,re_t,myData,logicalValidLFPs,LFPEndPts,LFPAmplitude,LFPAngles,fi,sourceElec,sourceReg,targetElecs,targetReg,well_spike_dyn,nYbin,thresh_mult,parent_dir,save_dir)

relationTable=table();
row=1;

% thetaAmpThresh=std(LFPAmplitude)*thresh_mult;
thetaAmpThresh=5;
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
for nElec=1:length(targetElecs)
    
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
    validWellBursts=[];
    for nBursts=1:size(myBurstBounds,1)
        validWellBursts=[validWellBursts,myBurstBounds(nBursts,1):myBurstBounds(nBursts,2)];
    end
    % well_spikes_idx=t(validWellSpikes);
    validWellBursts=round(remap(validWellBursts,1,length(t),1,length(re_t)));
    logicalValidBursts=zeros(1,length(re_t));
    logicalValidBursts(validWellBursts)=1;
    logicalValidBursts=logical(logicalValidBursts);

    nBurstsCounter(nElec)=sum(logicalBurstStarts & logicalValidLFPs);

    % figure
    wellBurstAngles=LFPAngles(logicalValidBursts & logicalValidLFPs);
    % wellBurstAngles=[wellBurstAngles-360,wellBurstAngles];
    wellBurstAmp=LFPAmplitude(logicalValidBursts & logicalValidLFPs);

    wellBurstStartAngles=LFPAngles(logicalBurstStarts & logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & logicalValidLFPs);

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_spb=well_spb(burstIdx);

    %% Amplitude MI
    fAmp=figure('Name',targetElecs(nElec)+" amp",'NumberTitle','off','WindowState','fullscreen');
    % tf=tiledlayout(3,2,"Padding","tight","TileSpacing","tight");
    % well spikes vs amplitude
    % nexttile
    histogram(wellBurstAmp,ampEdges)
    ampProbs=histcounts(wellBurstAmp,ampEdges,"Normalization","probability");
    ampCounts=histcounts(wellBurstAmp,ampEdges);
    ampMI=modulationIndex(ampProbs);
    ampPval=shuffleLFP_ModIdx(ampMI,LFPAmplitude,find(logicalValidLFPs),logicalValidBursts,ampEdges,nIter);
    % ampPval=circShuffleLFP_ModIdx(ampMI,LFPAmplitude,LFPEndPts,logicalValidSpikes,ampEdges,nIter);
    set(gca,"XScale","log")
    % set(gca,"FontSize",24)
    % axis square
    pbaspect([2,1,1])
    xlim([min(ampEdges),max(ampEdges)])
    xticks(ampEdges(1:4:end))
    xticklabels(round(ampEdges(1:4:end),2,"significant"))
    set(gca,'XMinorTick','off')
    ax=gca;
    ax.LineWidth=4;
    ax.TickLength=[0.01 0.05];

    title("Well Burst Times VS Amplitude")
    if ampPval<1/nIter
        subtitle("mod idx="+round(ampMI,2)+" p<"+1/nIter)
    else
        subtitle("mod idx="+round(ampMI,2,"significant")+" p="+ampPval)
    end
    ylabel("Burst Times")
    xlabel("Amplitude uV")
    set(gca,"FontSize",50)
    set(gca,"Position",[0.13,0.2,0.7750,0.6])
    
    % saveas(gcf,fullfile(save_dir,"AmpMI "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi+".png"),"png")
    %% Angle MI
    %well spikes vs angle
    % nexttile
    fAngle=figure('Name',targetElecs(nElec)+" angle",'NumberTitle','off','WindowState','fullscreen');
    histogram([wellBurstAngles-360,wellBurstAngles],angleEdges2Cycle)
    angleProbs=histcounts(wellBurstAngles,angleEdges,"Normalization","probability");
    angleCounts=histcounts(wellBurstAngles,angleEdges);
    angleMI=modulationIndex(angleProbs);
    anglePval=fftShuffleLFP_ModIdx(angleMI,myData,LFPEndPts,logicalValidBursts,angleEdges,nIter,"angle");
    % axis square
    pbaspect([2,1,1])
    xlim([-180,180])
    % xticks(angleEdges2Cycle(1:2:end))
    xticks([-180,-90,0,90,180])
    title("Well Burst Times VS Angle")
    if anglePval<1/nIter
        subtitle("mod idx="+round(angleMI,2)+" p<"+1/nIter)
    else
        subtitle("mod idx="+round(angleMI,2)+" p="+anglePval)
    end
    ylabel("BurstTimes")
    xlabel("Angle")
    xtickangle(0)
    set(gca,"FontSize",50)
    set(gca,"Position",[0.13,0.2,0.7750,0.6])
    ax=gca;
    ax.LineWidth=4;
    ax.TickLength=[0.01 0.05];

    hold on
    plot(linspace(-180,180,40),0.5*max(angleCounts)*sin(linspace(-pi,pi,40))+(0.5*max(angleCounts)),'r','LineWidth',4)
    hold off

    X=[[wellBurstAngles,wellBurstAngles-360];[wellBurstAmp,wellBurstAmp]]';
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

    else
        bincount_cells_xy{nElec}=[];
        binxcenters{nElec}=convert_edges_2_centers([-360:18:360]);
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1)));
        binxedges{nElec}=[-360:18:360];
        binyedges{nElec}=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);

        relationTable.fi(row)=fi;
        relationTable.sourceElec(row)=sourceElec;
        relationTable.sourceReg(row)=sourceReg;
        relationTable.targetElec(row)=targetElecs(nElec);
        relationTable.targetReg(row)=targetReg(nElec);
        relationTable.ampMI(row)=ampMI;
        relationTable.angleMI(row)=angleMI;
        relationTable.ampPval(row)=ampPval;
        relationTable.anglePval(row)=anglePval;
        relationTable.nAmpSpikesMax(row)=max(ampCounts);
        relationTable.nAngleSpikesMax(row)=max(angleCounts);
        relationTable.nHeatmapMax(row)=0;
        row=row+1;
        continue
    end
    % saveas(gcf,fullfile(save_dir,"AngleMI "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi),"png")

    relationTable.fi(row)=fi;
    relationTable.sourceElec(row)=sourceElec;
    relationTable.sourceReg(row)=sourceReg;
    relationTable.targetElec(row)=targetElecs(nElec);
    relationTable.targetReg(row)=targetReg(nElec);
    relationTable.ampMI(row)=ampMI;
    relationTable.angleMI(row)=angleMI;
    relationTable.ampPval(row)=ampPval;
    relationTable.anglePval(row)=anglePval;
    relationTable.nAmpSpikesMax(row)=max(ampCounts);
    relationTable.nAngleSpikesMax(row)=max(angleCounts);
    relationTable.nHeatmapMax(row)=max(bincount_cells_xy{nElec},[],"all");

    if max(ampCounts)>20 && max(angleCounts)>20 && max(bincount_cells_xy{nElec},[],"all")>6
        saveas(fAmp,fullfile(save_dir,"AmpMI "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi+".png"),"png")
        saveas(fAngle,fullfile(save_dir,"AngleMI "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi),"png")
    end

    row=row+1;
end

maxHeightAll=cellfun(@(x) max(x,[],"all"),bincount_cells_xy,'UniformOutput',false);
maxHeightAll(cellfun(@isempty,maxHeightAll))={0};
maxHeightAll=max(cell2mat(maxHeightAll));

row=1;

for nElec=1:length(targetElecs)
    if relationTable.nAmpSpikesMax(row)>20 && relationTable.nAngleSpikesMax(row)>20 && relationTable.nHeatmapMax(row)>6
        myFig=figure('Name',targetElecs(nElec)+" heatmap",'NumberTitle','off','WindowState','fullscreen');
        myFig.Position(3)=679.2;


        myCLim=[0,maxHeightAll];
        imagesc(binxcenters{nElec},binycenters{nElec},flipud(rot90(bincount_cells_xy{nElec},1)),myCLim)
        ylim([min(binyedges{nElec}),max(binyedges{nElec})])
        xlim([-180,180])

        xlabel("Angle (degrees)")
        ylabel("Amplitude uV")

        set(gca,"FontSize",50)

        xticks([-180,-90,0,90,180])
        yticks(binyedges{nElec}(1:4:end))
        yticklabels(round(binyedges{nElec}(1:4:end),2,"significant"))
        set(gca,'YDir','normal')
        set(gca,"YScale","log")

        cb = colorbar("FontSize",40);
        colormap hot
        cb.Limits=myCLim;
        clim(myCLim)
        ylabel(cb,'Burst ms','FontSize',40,'Rotation',270)
        if maxHeightAll>=5
            cb.Ticks=round(linspace(0,maxHeightAll,5));
        end
        % saveas(gcf,fullfile(save_dir,"HeatMap "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi),"png")
        axis square

        saveas(gcf,fullfile(save_dir,"HeatMap "+sourceElec+"-"+targetElecs(nElec)+" FID "+fi),"png")
    end
    row=row+1;
end

end