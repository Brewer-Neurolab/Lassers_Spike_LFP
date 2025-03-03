function miTable=underthresh_CMI_heatmap_formatter(jpxyz,jpxz,jpyz,bincounts_cells,binxcenters,binxedges,...
    binycenters,binyedges,binzcenters,binzedges,bincounts_x,bincounts_y,bincounts_z,...
    sourceElec,targetElecs,sourceProps,nBursts,nYbin)
% Formats colordata for heatmaps from a cell array where each cell contains
% a 2D matrix of color values. Flows and creates a tiled layout for cell
% array.

miTable=table();

maxHeightAll=cellfun(@(x) max(x,[],"all"),bincounts_cells,'UniformOutput',false);
maxHeightAll(cellfun(@isempty,maxHeightAll))={0};
maxHeightAll=max(cell2mat(maxHeightAll));

myCLim=[0,maxHeightAll];

figure
t=tiledlayout("flow","TileSpacing","tight","Padding","tight");
ax=[];

for nElec=1:length(bincounts_cells)
    % ax(nElec)=nexttile;
    if ~isempty(bincounts_cells{nElec}) %&& sum(bincounts_cells{nElec},"all")*2>=numel(binxcenters{nElec})*numel(binycenters{nElec})
        ax(nElec)=nexttile; %uncomment for only significant heatmaps
        %check significance of mod idx
        myMutInfo=condMutualInfo(jpxyz{nElec},jpxz{nElec},jpyz{nElec},...
            bincounts_x{nElec}/sum(bincounts_x{nElec}),...
            bincounts_y{nElec}/sum(bincounts_y{nElec}),...
            bincounts_z{nElec}/sum(bincounts_z{nElec}));
        if myMutInfo<=0.01
            myMutInfo=0;
        end
        pval=underthresh_shuffleCMI(myMutInfo,sourceProps.LFPAmps,sourceProps.LFPAngles,sourceProps.logicalValidLFPs,...
            sourceProps.well_spike_dyn,targetElecs(nElec),sourceProps.fi,sourceProps.t,sourceProps.re_t,sourceProps.nIter,nYbin);       
        imagesc(binxcenters{nElec},binycenters{nElec},flipud(rot90(bincounts_cells{nElec},1)),myCLim)
        ylim([min(binyedges{nElec}),max(binyedges{nElec})])
        if pval<1/sourceProps.nIter
            title(targetElecs(nElec)+" CMI="+round(myMutInfo,2,"significant")+" p<"+1/sourceProps.nIter+" n="+nBursts(nElec),"FontSize",12)
        else
            title(targetElecs(nElec)+" CMI="+round(myMutInfo,2,"significant")+" p="+pval+" n="+nBursts(nElec),"FontSize",12)
        end
        set(gca,"FontSize",14)
        xticks(binxedges{nElec})
        yticks(binyedges{nElec}(1:2:end))
        yticklabels(round(binyedges{nElec}(1:2:end),2,"significant"))
        set(gca,'YDir','normal')
        set(gca,"YScale","log")
        axis("square")
    else
        myMutInfo="NA";
        pval="NA";
        % title(targetElecs(nElec)+" MI="+myMutInfo+" p="+pval+" n="+nBursts(nElec),"FontSize",12)
    end

    miTable.fi(nElec)=sourceProps.fi;
    miTable.sourceElec(nElec)=sourceElec;
    miTable.targetElec(nElec)=targetElecs(nElec);
    miTable.MI(nElec)=myMutInfo;
    miTable.pval(nElec)=pval;
    disp(sourceElec+" "+nElec+" of "+length(bincounts_cells)+" MI done.")
end

cb = colorbar("FontSize",18);
colormap hot
cb.Layout.Tile = 'east';
cb.Limits=myCLim;
clim(myCLim)
cb.Ticks=round(linspace(0,maxHeightAll,10));
% set(cb,"ColorScale","log")

end