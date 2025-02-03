function scatter_BL_v_SPB(well_spike_dyn, MI_tbl)

figure
t=tiledlayout("flow","TileSpacing","tight","Padding","tight");
ax=[];

for nElec=1:height(MI_tbl)
    ax(nElec)=nexttile;

    burst_length=well_spike_dyn.BurstDuration{well_spike_dyn.fi==MI_tbl.fi(nElec)&well_spike_dyn.channel_name==MI_tbl.targetElec(nElec)};
    spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==MI_tbl.fi(nElec)&well_spike_dyn.channel_name==MI_tbl.targetElec(nElec)};

    scatter(burst_length,spb)

    
    xlabel("Log Burst Duration (ms)")
    ylabel("Spikes Per Burst")

    hold on 
    % lsline
    mdl=fitlm(burst_length,spb);
    x=linspace(min(xlim),max(xlim),10);
    y=x.*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1);
    plot(x,y)
    hold off

    title(MI_tbl.targetElec(nElec)+" r^2="+round(mdl.Rsquared.Adjusted,2))

    set(gca,"FontSize",16)
end

end