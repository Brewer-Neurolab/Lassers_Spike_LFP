function [pval,MI_Vec]=underthresh_shuffleCMI(myMI,LFPAmps,LFPAngles,logicalValidLFPs,well_spike_dyn,myElec,fi,t,re_t,nIter,nYbin)
MI_Vec=[];
pval=[];

well_burst_bounds=well_spike_dyn.BurstBounds{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==myElec};
well_burst_starts=well_burst_bounds(:,1);
% remap burst starts to new sampling rate
well_burst_starts=round(remap(well_burst_starts,1,length(t),1,length(re_t)));
logicalBurstStarts=zeros(1,length(re_t));
logicalBurstStarts(well_burst_starts)=1;
thetaAmpThresh=std(LFPAmps);

rng('default')
nPermutations=nIter;
shuffledLFPAngles=zeros(nPermutations,length(LFPAngles));
shuffledLFPAmps=zeros(nPermutations,length(LFPAmps));

%calculate MI for each permuatation
for nPerm=1:nPermutations
    randPermAngles=randperm(length(LFPAngles));
    shuffledLFPAngles(nPerm,:)=LFPAngles(randPermAngles);
    randPermAmps=randperm(length(LFPAmps));
    shuffledLFPAmps(nPerm,:)=LFPAmps(randPermAmps);

    wellBurstStartAngles=shuffledLFPAngles(nPerm,logicalBurstStarts & ~logicalValidLFPs);
    wellBurstStartAmp=shuffledLFPAmps(nPerm,logicalBurstStarts & ~logicalValidLFPs);

    %repeat for spikes per burst
    repwellBurstStartAngles=[];
    repwellBurstStartAmp=[];
    repwellSPB=[];

    burstIdx=ismembertol(well_burst_starts,find(logicalBurstStarts & ~logicalValidLFPs),1e-10);
    burstIdx=find(burstIdx);
    well_spb=well_spike_dyn.SpikeperBurst{well_spike_dyn.fi==fi & well_spike_dyn.channel_name==myElec};

    for nBursts=1:length(wellBurstStartAngles)
        % repwellBurstStartAngles=[repwellBurstStartAngles,repmat(wellBurstStartAngles(nBursts),1,well_spb(burstIdx(nBursts)))];
        % repwellBurstStartAmp=[repwellBurstStartAmp,repmat(wellBurstStartAmp(nBursts),1,well_spb(burstIdx(nBursts)))];
        % repwellSPB=[repwellSPB,repmat(well_spb(burstIdx(nBursts)),1,well_spb(burstIdx(nBursts)))];
        repwellBurstStartAngles=[repwellBurstStartAngles,wellBurstStartAngles(nBursts)];
        repwellBurstStartAmp=[repwellBurstStartAmp,wellBurstStartAmp(nBursts)];
        repwellSPB=[repwellSPB,well_spb(burstIdx(nBursts))];
    end

    nbins_x=8;
    nbins_z=40;

    X=[repwellBurstStartAngles;repwellBurstStartAmp]';
    edges={linspace(0,360,nbins_x+1),logspace(log10(min(LFPAmps)),log10(thetaAmpThresh),nYbin+1)};

    if ~isempty(repwellBurstStartAngles)

        %calculates pxy
        [N]=hist3(X,'Edges',edges,'CDataMode','manual','FaceColor','interp');
        bincount_cells_xy=N(1:nbins_x,1:nYbin);
        binxcenters=convert_edges_2_centers(linspace(0,360,nbins_x+1));
        binycenters=10.^convert_edges_2_centers(log10(logspace(log10(min(LFPAmps)),log10(thetaAmpThresh),nYbin+1)));
        binzcenters=convert_edges_2_centers(linspace(0,max(repwellSPB),nbins_z));

        binxedges=linspace(0,360,nbins_x+1);
        binyedges=logspace(log10(min(repwellBurstStartAmp)),log10(max(LFPAmps)),nYbin+1);
        binzedges=linspace(0,max(repwellSPB),nbins_z);

        %calculate px
        bincount_cells_x=histcounts(repwellBurstStartAngles,linspace(0,360,nbins_x+1));

        %calculate py
        bincount_cells_y=histcounts(repwellBurstStartAmp,logspace(log10(min(LFPAmps)),log10(thetaAmpThresh),nYbin+1));

        %calculate pz
        bincount_cells_z=histcounts(repwellSPB,binzedges);

        [~,~,bins1]=histcounts(repwellBurstStartAngles,binxedges);
        [~,~,bins2]=histcounts(repwellBurstStartAmp,binyedges);
        [~,~,bins3]=histcounts(repwellSPB,binzedges);

        jpxyz=zeros([nbins_x,nYbin,nbins_z]);
        for x=1:nbins_x
            for y=1:nYbin
                for z=1:nbins_z
                    jpxyz(x,y,z)=sum(bins1==x & bins2==y & bins3==z);
                end
            end
        end
        jpxyz=jpxyz./sum(jpxyz,"all");

        jpxz=zeros([nbins_x,nbins_z]);
        for x=1:nbins_x
            for z=1:nbins_z
                jpxz(x,z)=sum(bins1==x & bins3==z);
            end
        end
        jpxz=jpxz./sum(jpxz,"all");

        jpyz=zeros([nYbin,nbins_z]);
        for y=1:nYbin
            for z=1:nbins_z
                jpyz(y,z)=sum(bins2==y & bins3==z);
            end
        end
        jpyz=jpyz./sum(jpyz,"all");

    else
        binxcenters=[];
        binycenters=[];
        binzcenters=[];
        binxedges=[];
        binyedges=[];
        binzedges=[];
        jpxyz=[];
        jpxz=[];
        jpyz=[];
    end

    MI_Vec(nPerm)=condMutualInfo(jpxyz,jpxz,jpyz,...
            bincount_cells_x/sum(bincount_cells_x),...
            bincount_cells_y/sum(bincount_cells_y),...
            bincount_cells_z/sum(bincount_cells_z));

    %comment out when not testing
    % figure
    % imagesc(binxcenters,binycenters,rot90(flipud(bincounts_cells),-1))
    % xticks(binxcenters)
    % yticks(binycenters)
    % yticklabels(round(binycenters,2,"significant"))
    % set(gca,'YDir','normal')
    % set(gca,"YScale","log")
    % axis("square")
    % cb = colorbar;
    % cb.Limits=[0,200];
    % colormap hot
    % clim([0,200])
    % 
end

%calculate total pval
if ~isempty(MI_Vec)
    pval=sum(MI_Vec>=myMI)/length(MI_Vec);
end

% for testing and validating pval
% figure
% histogram(MI_Vec)
% hold on
% xline(myMI)
% hold off

end