%% Axon to Well Amplitude and Spike Relation Test for Gamma
%% Setup
clear
clc
close all

saveDir="C:\BrewerLabResearch\GitHub\Lassers_Spike_LFP\Gamma Scripts";

parent_axons_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\downsampled tunnels\Low_Gamma"; %change when switching
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

parent_wells_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes 5SD Min";
wells_dir=dir(parent_wells_dir);
wells_folders=string({wells_dir.name});
wells_folders=wells_folders([wells_dir.isdir]);
wells_folders=wells_folders(3:end);

axon_spikes=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axon_spikes=axon_spikes.allregion_unit_matched_stim;

interRegions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];
subregions=["EC","DG","CA3","CA1"];
ff_axon_tbl=table();

row=1;
for fi=1:length(axon_spikes)
    for nelec=1:height(axon_spikes{fi})
        if ~isempty(axon_spikes{fi}.up_ff{nelec})
            ff_axon_tbl.fi(row)=fi;
            ff_axon_tbl.Subregion(row)=axon_spikes{fi}.Subregion(nelec);
            ff_axon_tbl.interRegi(row)=find(interRegions==axon_spikes{fi}.Subregion(nelec));
            electrodes=axon_spikes{fi}.("Electrode Pairs")(nelec);
            electrodes=split(electrodes,{'-'});
            ff_axon_tbl.Electrode(row)=electrodes(1);
            tunnelReg=split(axon_spikes{fi}.Subregion(nelec),{'-'});
            ff_axon_tbl.FFReg(row)=tunnelReg(2);
            ff_axon_tbl.subi(row)=find(subregions==tunnelReg(2));
            row=row+1;
        end
    end
end

%3.5 SD
%%well_spike_dyn=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat")

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;
%% Compute MI and heat maps
clc

%testing range
testing_idx=[find(ff_axon_tbl.Subregion=="CA3-CA1")]';

relationTable=[];

for nFF=1:height(ff_axon_tbl)
    data=load(fullfile(parent_axons_dir,axons_folders(ff_axon_tbl.fi(nFF)),ff_axon_tbl.Electrode(nFF)+".mat"));
    % re_fs=data.re_fs;
    re_fs=1000;
    fs=25000;
    data=data.filtered_data;
    t_rec=300;
    re_t=0:1/re_fs:t_rec-(1/re_fs);
    t=0:1/fs:t_rec-(1/fs);

    % plot axon data tagged
    thresh_mult=1;

    %define max number of samples for combining LFPs
    nsamples_combine_thresh=(1/100)*re_fs*3; %change when switching
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest gamma cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/100)*minLFPCycles*re_fs; %change when switching

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    validLFPIndex=[];
    for nEndPts=1:size(LFPEndPts,1)
        validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
    end
    logicalValidLFPs=ones(1,length(re_t));
    %%logicalValidLFPs(validLFPIndex)=1;

    % Regression tests
    targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF));
    targetReg=subregions(well_spike_dyn.regi(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF)));
    sourceReg=ff_axon_tbl.Subregion(nFF);
    
    myTable=sourceLFP_targetSpike_relations_NoThresh(t,re_t,data,logicalValidLFPs,LFPEndPts,LFPAmplitude,LFPAngles,ff_axon_tbl.fi(nFF),ff_axon_tbl.Electrode(nFF),sourceReg,targetElecs,targetReg,well_spike_dyn,20,thresh_mult,...
        fullfile(parent_wells_dir,wells_folders(ff_axon_tbl.fi(nFF))+"\"),...
        "C:\BrewerLabResearch\GitHub\Lassers_Spike_LFP\Gamma Scripts\Images_Low"); %change when switching


    if isempty(relationTable)
        relationTable=myTable;
    else
        relationTable=[relationTable;myTable];
    end

    disp(nFF+" of "+height(ff_axon_tbl))

    close all force
end

save(fullfile(saveDir,"relationTable"),"relationTable")
%% Good Relationships

%calculate false discovery rates for amp and angle
[~,IAmp]=sort(relationTable.ampPval);
[relationTable.ampFDR,relationTable.ampFDRh,largestAmpP]=myFalseDiscoveryRate(IAmp,length(IAmp),0.05,relationTable.ampPval);

[~,IAngle]=sort(relationTable.anglePval);
[relationTable.angleFDR,relationTable.angleFDRh,largestAngleP]=myFalseDiscoveryRate(IAngle,length(IAngle),0.05,relationTable.anglePval);

goodRelationsTbl=relationTable((relationTable.ampPval<0.05 | relationTable.anglePval<0.05)...
    & relationTable.nAmpSpikesMax>20 & relationTable.nAngleSpikesMax>20 & relationTable.nHeatmapMax>10,:);
goodRelationsTblFDR=relationTable((relationTable.ampFDRh & relationTable.angleFDRh)...
    & relationTable.nAmpSpikesMax>20 & relationTable.nAngleSpikesMax>20 & relationTable.nHeatmapMax>10,:);

save(fullfile(saveDir,"goodRelationsTbl"),"goodRelationsTbl")
save(fullfile(saveDir,"goodRelationsTblFDR"),"goodRelationsTblFDR")

