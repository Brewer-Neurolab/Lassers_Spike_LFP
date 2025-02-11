function scatter_BL_v_SPB(well_spike_dyn, MI_tbl,t,re_t,logicalValidLFPs)

figure
tiledlayout("flow","TileSpacing","tight","Padding","tight");
ax=[];

for nElec=1:height(MI_tbl)
    ax(nElec)=nexttile;

    burst_length=well_spike_dyn.BurstDuration{well_spike_dyn.fi==MI_tbl.fi(nElec)&well_spike_dyn.channel_name==MI_tbl.targetElec(nElec)};
    spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==MI_tbl.fi(nElec)&well_spike_dyn.channel_name==MI_tbl.targetElec(nElec)};

    scatter(burst_length,spb,".b")

    
    xlabel("Burst Duration (ms)")
    ylabel("Spikes Per Burst")

    hold on 
    % lsline
    mdl1=fitlm(burst_length,spb);
    x=linspace(min(xlim),max(xlim),10);
    y=x.*mdl1.Coefficients.Estimate(2)+mdl1.Coefficients.Estimate(1);
    plot(x,y,'--b')
    hold off

    title(MI_tbl.targetElec(nElec))

    set(gca,"FontSize",16)

    % find ok bursts
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==MI_tbl.fi(nElec) & well_spike_dyn.channel_name==MI_tbl.targetElec(nElec)};
    well_burst_starts=well_burst_bounds(:,1);
    % remap burst starts to new sampling rate
    well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
    logicalBurstStarts=zeros(1,length(re_t));
    logicalBurstStarts(well_burst_starts)=1;

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);

    hold on
    % lsline
    mdl2=fitlm(burst_length(burstIdx),spb(burstIdx));
    x=linspace(min(xlim),max(xlim),10);
    y=x.*mdl2.Coefficients.Estimate(2)+mdl2.Coefficients.Estimate(1);
    plot(x,y,'--r')

    scatter(burst_length(burstIdx),spb(burstIdx),".r")
    hold off

    lgd=legend('',"All r^2="+round(mdl1.Rsquared.Adjusted,2)+newline+"Mean BL="+newline+mean(burst_length)+newline+"Mean SBP="+mean(spb)+newline+"SD BL="+std(burst_length)+newline+"SD SPB="+std(spb),...
        "High Amp r^2="+round(mdl2.Rsquared.Adjusted,2)+newline+"Mean BL="+mean(burst_length(burstIdx))+newline+"Mean SBP="+mean(spb(burstIdx))+newline+"SD BL="+std(burst_length(burstIdx))+newline+"SD SPB="+std(spb(burstIdx)),'');

    lgd.Location="eastoutside";

    axis square
end

end