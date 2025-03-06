function miTable=underthresh_CMI_axon_well_burst_start(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,fi,sourceElec,targetElecs,well_spike_dyn,nYbin,thresh_mult)

bincount_cells_x=[];
bincount_cells_y=[];
bincount_cells_z=[];
binxcenters=[];
binxedges=[];
binycenters=[];
binyedges=[];
binzcenters=[];
binzedges=[];

jpxyz=[];
jpxz=[];
jpyz=[];

for nElec=1:length(targetElecs)
    disp("Calculating "+targetElecs(nElec))
    well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};
    well_burst_starts=well_burst_bounds(:,1);
    % remap burst starts to new sampling rate
    well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
    logicalBurstStarts=zeros(1,length(re_t));
    logicalBurstStarts(well_burst_starts)=1;

    nBurstsCounter(nElec)=sum(logicalBurstStarts & ~logicalValidLFPs);

    % figure
    wellBurstStartAngles=LFPAngles(logicalBurstStarts & ~logicalValidLFPs);
    % wellBurstStartAngles=[wellBurstStartAngles-360,wellBurstStartAngles];
    wellBurstStartAmp=LFPAmplitude(logicalBurstStarts & ~logicalValidLFPs);

    thetaAmpThresh=std(LFPAmplitude)*thresh_mult;

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];
    repwellSPB=[];

    countwellBurstStartAngles=[];
    countwellBurstStartAmp=[];
    countwellSPB=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & ~logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==targetElecs(nElec)};

    for nBursts=1:length(wellBurstStartAngles)
        repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
        repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
        repwellSPB=[repwellSPB,repmat(well_spb(burstIdx(nBursts)),1,well_spb(burstIdx(nBursts)))];
        countwellBurstStartAngles=[countwellBurstStartAngles,wellBurstStartAngles(nBursts)];
        countwellBurstStartAmp=[countwellBurstStartAmp,wellBurstStartAmp(nBursts)];
        countwellSPB=[countwellSPB,well_spb(burstIdx(nBursts))];
    end

    nbins_x=8;
    nbins_z=40;

    X=[repwellBurstStartAngles;repwellBurstStartAmp]';
    edges={linspace(0,360,nbins_x+1),logspace(log10(min(LFPAmplitude)),log10(thetaAmpThresh),nYbin+1)};

    bincount_cells_xy=[];

    if ~isempty(repwellBurstStartAngles)

        %calculates pxy
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells_xy{nElec}=N(1:nbins_x,1:nYbin);
        binxcenters{nElec}=convert_edges_2_centers(linspace(0,360,nbins_x+1));
        binycenters{nElec}=10.^convert_edges_2_centers(log10(logspace(log10(min(LFPAmplitude)),log10(thetaAmpThresh),nYbin+1)));
        binzcenters{nElec}=convert_edges_2_centers(linspace(0,max(countwellSPB),nbins_z));

        binxedges{nElec}=linspace(0,360,nbins_x+1);
        binyedges{nElec}=logspace(log10(min(LFPAmplitude)),log10(thetaAmpThresh),nYbin+1);
        binzedges{nElec}=linspace(0,max(countwellSPB),nbins_z);

        %calculate px
        bincount_cells_x{nElec}=histcounts(countwellBurstStartAngles,linspace(0,360,nbins_x+1));

        %calculate py
        bincount_cells_y{nElec}=histcounts(countwellBurstStartAmp,logspace(log10(min(LFPAmplitude)),log10(thetaAmpThresh),nYbin+1));

        %calculate pz
        bincount_cells_z{nElec}=histcounts(countwellSPB,binzedges{nElec});

        [~,~,bins1]=histcounts(countwellBurstStartAngles,binxedges{nElec});
        [~,~,bins2]=histcounts(countwellBurstStartAmp,binyedges{nElec});
        [~,~,bins3]=histcounts(countwellSPB,binzedges{nElec});

        jpxyz{nElec}=zeros([nbins_x,nYbin,nbins_z]);
        for x=1:nbins_x
            for y=1:nYbin
                for z=1:nbins_z
                    jpxyz{nElec}(x,y,z)=sum(bins1==x & bins2==y & bins3==z);
                end
            end
        end
        jpxyz{nElec}=jpxyz{nElec}./sum(jpxyz{nElec},"all");

        jpxz{nElec}=zeros([nbins_x,nbins_z]);
        for x=1:nbins_x
            for z=1:nbins_z
                jpxz{nElec}(x,z)=sum(bins1==x & bins3==z);
            end
        end
        jpxz{nElec}=jpxz{nElec}./sum(jpxz{nElec},"all");

        jpyz{nElec}=zeros([nYbin,nbins_z]);
        for y=1:nYbin
            for z=1:nbins_z
                jpyz{nElec}(y,z)=sum(bins2==y & bins3==z);
            end
        end
        jpyz{nElec}=jpyz{nElec}./sum(jpyz{nElec},"all");

    else
        binxcenters{nElec}=[];
        binycenters{nElec}=[];
        binzcenters{nElec}=[];
        binxedges{nElec}=[];
        binyedges{nElec}=[];
        binzedges{nElec}=[];
        jpxyz{nElec}=[];
        jpxz{nElec}=[];
        jpyz{nElec}=[];
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
miTable=underthresh_CMI_heatmap_formatter(jpxyz,jpxz,jpyz,bincount_cells_xy,binxcenters,binxedges,...
    binycenters,binyedges,binzcenters,binzedges,bincount_cells_x,bincount_cells_y,...
    bincount_cells_z,sourceElec,targetElecs,sourceProps,nBurstsCounter,nYbin);

end