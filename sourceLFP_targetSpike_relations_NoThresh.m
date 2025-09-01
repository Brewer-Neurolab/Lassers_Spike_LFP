function relationTable=sourceLFP_targetSpike_relations_NoThresh(t,re_t,myData,logicalValidLFPs,LFPEndPts,LFPAmplitude,LFPAngles,fi,sourceElec,sourceReg,targetElecs,targetReg,well_spike_dyn,nYbin,thresh_mult,parent_dir,save_dir)

logicalValidLFPs=logical(logicalValidLFPs);
relationTable=table();
row=1;

thetaAmpThresh=5;
ampEdges=logspace(log10(thetaAmpThresh),log10(max(LFPAmplitude)),nYbin+1);
% ampEdges=logspace(log10(min(LFPAmplitude(logicalValidLFPs))),log10(max(LFPAmplitude)),nYbin+1);
ampCenters=convert_edges_2_centers(log10(ampEdges));
ampCenters=10.^ampCenters;
angleEdges=[0:18:360];
angleCenters=convert_edges_2_centers(angleEdges);
angleEdges2Cycle=[-360:18:360];

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

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_spb=well_spb(burstIdx);

    %% Amplitude MI
    fAmp=figure('Name',targetElecs(nElec)+" amp",'NumberTitle','off','WindowState','fullscreen');
    fAmp_axis=axes('Parent',fAmp);
    % tf=tiledlayout(3,2,"Padding","tight","TileSpacing","tight");
    % well spikes vs amplitude
    % nexttile
    histogram(wellBurstAmp,ampEdges)
    ampProbs=histcounts(wellBurstAmp,ampEdges,"Normalization","probability");
    ampCounts=histcounts(wellBurstAmp,ampEdges);
    ampMI=modulationIndex(ampProbs);
    
    %shuffling tests for pval
    % ampPval=shuffleLFP_ModIdx(ampMI,LFPAmplitude,find(logicalValidLFPs),logicalValidSpikes,ampEdges,nIter);
    % ampPval=circShuffleLFP_ModIdx(ampMI,LFPAmplitude,LFPEndPts,logicalValidSpikes,ampEdges,nIter);
    % ampPval=jitterLFP_ModIdx(ampMI,LFPAmplitude,find(logicalValidLFPs),logicalValidSpikes,ampEdges,nIter,100,1000);

    %slope fit test for pval
    % mdl=fitlm(log10(ampCenters),hamp.Values);
    % ampPval=mdl.Coefficients.pValue(2);
    % [coeff,stats,optlim,~,~]=find_semipowerlawfit_with_grid_search(ampCenters,ampCounts,[min(ampCenters),max(ampCenters)]);
    % ampPval=stats.pValue(2);
    % rsq=stats.Rsquared;
    % if ~isempty(optlim)
    %     x=linspace(optlim(1),optlim(2),20);
    %     y=coeff(2).*log10(x)+coeff(1);
    %     hold(fAmp_axis,"on")
    %     plot(fAmp_axis,x,y,'r--')
    %     hold(fAmp_axis,"off")
    % end
    %
    % if isempty(coeff)
    %     coeff=[NaN,NaN];
    % end
    if max(ampCounts)>20
        [coeff,stats,optlim,~,~]=find_powerlawfit_with_grid_search(ampCenters,ampCounts,[min(ampCenters),max(ampCenters)]);
        ampPval=stats.pValue(2);
        rsq=stats.Rsquared;
        if ~isempty(optlim)
            x=linspace(optlim(1),optlim(2),20);
            y=10.^(coeff(2).*log10(x)+coeff(1));
            hold(fAmp_axis,"on")
            plot(fAmp_axis,x,y,'r--','LineWidth',4)
            hold(fAmp_axis,"off")
            xlim([min(ampEdges),max(ampEdges)])
            % xticks(ampEdges(1:4:end))
            % xticks(logspace(0,5,12))
            % xticklabels(round(ampEdges(1:4:end),2,"significant"))
            % xticklabels(logspace(0,5,12))
            xticks([5,10,20,50,100,200,500,1000,5000])
            xticklabels([5,10,20,50,100,200,500,1000,5000])
            set(fAmp_axis,'XMinorTick','on')

            yticks([2,5,10,20,50,100,200,500,1000,5000])
            yticklabels([2,5,10,20,50,100,200,500,1000,5000])
            set(fAmp_axis,'YMinorTick','on')
            % yticks(logspace(0,max(ampCounts),12))
            % yticklabels(round(logspace(0,max(ampCounts),12),2,"significant"))
        end

        if isempty(coeff)
            coeff=[NaN,NaN];
        end
    else
        coeff=[NaN,NaN];
        rsq=NaN;
        ampPval=NaN;
    end
    set(fAmp_axis,"XScale","log")
    set(fAmp_axis,"YScale","log")
    % set(gca,"FontSize",24)
    % axis square
    pbaspect([2,1,1])
    

    ax=fAmp_axis;
    ax.LineWidth=4;
    ax.TickLength=[0.01 0.05];

    title("Well Spikes VS Amplitude")
    % if ampPval<1/nIter
    %     subtitle("R^2="+round(rsq,2)+" p<"+1/nIter)
    % else
    %     subtitle("R^2="+round(rsq,2,"significant")+" Slope "+coeff(2)+" p="+ampPval)
    % end

    subtitle("R^2="+round(rsq,2,"significant")+" Slope "+round(coeff(2),2,'significant')+" p="+round(ampPval,2,"significant"))

    ylabel("Spikes")
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
    anglePval=fftShuffleLFP_ModIdx(angleMI,myData,LFPEndPts,logicalValidSpikes,angleEdges,nIter,"angle");
    % axis square
    pbaspect([2,1,1])
    xlim([-180,180])
    % xticks(angleEdges2Cycle(1:2:end))
    xticks([-180,-90,0,90,180])
    title("Well Spikes VS Angle")
    if anglePval<1/nIter
        subtitle("mod idx="+round(angleMI,2)+" p<"+1/nIter)
    else
        subtitle("mod idx="+round(angleMI,2)+" p="+anglePval)
    end
    ylabel("Spikes")
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
        relationTable.rsq(row)=rsq;
        relationTable.slope(row)=coeff(2);
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
    relationTable.rsq(row)=rsq;
    relationTable.slope(row)=coeff(2);
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
        yticks([2,5,10,20,50,100,200,500,1000,5000])
        yticklabels([2,5,10,20,50,100,200,500,1000,5000])
        set(gca,'YDir','normal')
        set(gca,"YScale","log")

        cb = colorbar("FontSize",40);
        colormap hot
        cb.Limits=myCLim;
        clim(myCLim)
        ylabel(cb,'#Spikes','FontSize',40,'Rotation',270)
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