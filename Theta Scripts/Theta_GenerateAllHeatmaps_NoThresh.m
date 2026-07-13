%% Axon to Well Amplitude and Spike Relation Test for Theta
%% Setup
clear
clc
close all

saveDir="C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Theta Scripts";

parent_axons_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels\theta";
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

parent_wells_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min";
wells_dir=dir(parent_wells_dir);
wells_folders=string({wells_dir.name});
wells_folders=wells_folders([wells_dir.isdir]);
wells_folders=wells_folders(3:end);

axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axon_spikes=axon_spikes.allregion_unit_matched_stim;

interRegions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];
subregions=["EC","DG","CA3","CA1"];
ff_axon_tbl=table();

row=1;
for fi=1:length(axon_spikes)
    for nelec=1:height(axon_spikes{fi})
        if ~isempty(axon_spikes{fi}.up_ff{nelec}) | (~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC")
            ff_axon_tbl.fi(row)=fi;
            ff_axon_tbl.Subregion(row)=axon_spikes{fi}.Subregion(nelec);
            ff_axon_tbl.interRegi(row)=find(interRegions==axon_spikes{fi}.Subregion(nelec));
            electrodes=axon_spikes{fi}.("Electrode Pairs")(nelec);
            electrodes=split(electrodes,{'-'});
            ff_axon_tbl.Electrode(row)=electrodes(1);
            tunnelReg=split(axon_spikes{fi}.Subregion(nelec),{'-'});
            ff_axon_tbl.FFReg(row)=tunnelReg(2);
            ff_axon_tbl.subi(row)=find(subregions==tunnelReg(2));

            if isempty(axon_spikes{fi}.up_ff{nelec}) & ~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC"
                ff_axon_tbl.Subregion(row)="EC-CA1";
                ff_axon_tbl.FFReg(row)="CA1";
                ff_axon_tbl.subi(row)=4;
            end

            row=row+1;
        end
    end
end

%3.5 SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

all_reg=[interRegions,"EC-CA1"];

%% Compute MI and heat maps
clc

%testing range
testing_idx=[find(ff_axon_tbl.Subregion=="CA3-CA1")]';

relationTable=[];

for nFF=27%1:height(ff_axon_tbl)
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
    nsamples_combine_thresh=(1/10)*re_fs*3;
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest theta cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/10)*minLFPCycles*re_fs;

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    validLFPIndex=[];
    for nEndPts=1:size(LFPEndPts,1)
        validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
    end
    logicalValidLFPs=zeros(1,length(re_t)); % uncomment for lower bound LFP
    logicalValidLFPs=ones(1,length(re_t)); % considers all LFPs
    % logicalValidLFPs(validLFPIndex)=1;

    % Regression tests
    targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF));
    targetReg=subregions(well_spike_dyn.regi(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF)));
    sourceReg=ff_axon_tbl.Subregion(nFF);
    myTable=sourceLFP_targetSpike_relations_NoThresh(t,re_t,data,logicalValidLFPs,LFPEndPts,LFPAmplitude,LFPAngles,ff_axon_tbl.fi(nFF),ff_axon_tbl.Electrode(nFF),sourceReg,targetElecs,targetReg,well_spike_dyn,20,thresh_mult,...
        fullfile(parent_wells_dir,wells_folders(ff_axon_tbl.fi(nFF))+"\"),...
        "C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Images\Theta\Spikes No Thresh LogLog 5SD 1000 iter");

    if isempty(relationTable)
        relationTable=myTable;
    else
        relationTable=[relationTable;myTable];
    end

    disp(nFF+" of "+height(ff_axon_tbl))

    close all force
end

save(fullfile(saveDir,"relationTable_NoThresh"),"relationTable")
%% Good Relationships

load(fullfile(saveDir,"relationTable_NoThresh"),"relationTable")

goodRelationsTbl=table();
goodRelationsTblFDR=table();
largestAngleP=[];

for i=1:length(all_reg)
    %calculate false discovery rates for amp and angle
    sub_tbl=relationTable(relationTable.sourceReg==all_reg(i),:);

    % [~,IAmp]=sort(relationTable.ampPval);
    % [sub_tbl.ampFDR,sub_tbl.ampFDRh,largestAmpP]=myFalseDiscoveryRate(IAmp,length(IAmp),0.05,sub_tbl.ampPval);

    [~,IAngle]=sort(sub_tbl.anglePval);
    [sub_tbl.angleFDR,sub_tbl.angleFDRh,largestAngleP(i)]=myFalseDiscoveryRate(IAngle,length(IAngle),0.05,sub_tbl.anglePval);

    temp_tab1=sub_tbl((sub_tbl.ampPval<0.05 | sub_tbl.anglePval<0.05)...
        & sub_tbl.nAmpSpikesMax>20 & sub_tbl.nAngleSpikesMax>20 & sub_tbl.nHeatmapMax>10,:);
    goodRelationsTbl=[goodRelationsTbl;temp_tab1];

    % temp_tab2=sub_tbl((sub_tbl.angleFDRh)...
    %     & sub_tbl.nAmpSpikesMax>20 & sub_tbl.nAngleSpikesMax>20 & sub_tbl.nHeatmapMax>10,:);
    temp_tab2=sub_tbl((sub_tbl.angleFDRh& sub_tbl.nAmpSpikesMax>20 & sub_tbl.nAngleSpikesMax>20),:);
    goodRelationsTblFDR=[goodRelationsTblFDR;temp_tab2];
end

save(fullfile(saveDir,"goodRelationsTbl_NoThresh"),"goodRelationsTbl")
save(fullfile(saveDir,"goodRelationsTblFDR_NoThresh"),"goodRelationsTblFDR")

%% Avg Spike Angle By subregion
close all

load(fullfile(saveDir,"goodRelationsTbl_NoThresh"))
load(fullfile(saveDir,"goodRelationsTblFDR_NoThresh"))

% posGoodRel=goodRelationsTblFDR(goodRelationsTblFDR.slope>0,:);
posGoodRel=goodRelationsTblFDR;

% all_reg=[interRegions,"EC-CA1"];

sumAngleCounts=[];

for i=1:length(all_reg)
    figure
    angleProbs=cell2mat(posGoodRel.angleProbs(posGoodRel.sourceReg==all_reg(i)));
    angleCounts=cell2mat(posGoodRel.angleCounts(posGoodRel.sourceReg==all_reg(i)));
    nPairs(i)=size(angleProbs,1);
    meanAngleProbs=mean(angleProbs,1);
    stdAngleProbs=std(angleProbs,[],1);
    sumAngleCounts{i}=sum(angleCounts,1);
    histogram("BinEdges",[-180:18:180],"BinCounts",meanAngleProbs)
    hold on
    errorbar(convert_edges_2_centers([-180:18:180]),meanAngleProbs,stdAngleProbs,"LineStyle","none","Color","k")
    % title(all_reg(i))
    xlabel("Angle")
    ylabel("Spike Probability")
    % ylim([0,max(meanAngleProbs+stdAngleProbs)*1.1])
    ylim([0,0.2])
    xlim([-180,180])
    xticks(-180:36:180)
    set(gca,"FontSize",18)
    hold on
    x=-180:180;
    plot(-180:180,(sind(x)*0.1)+0.1,'r','LineWidth',2)
    originalMI(i)=modulationIndex(meanAngleProbs);

    niter=1000;
    shuffleMI=[];
    m=[];
    for j=1:niter
        permMat=[];
        for k=1:size(angleProbs,1)
            permMat(k,:)=randperm(size(angleProbs,2));
        end
        m(j,:)=mean(angleProbs(permMat),1)/sum(mean(angleProbs(permMat),1));
        shuffleMI(j)=modulationIndex(m(j,:));
    end

    pval(i)=sum(shuffleMI>originalMI(i))/niter;
end

%% ANOVA of mean spike angle
bin_centers=deg2rad(convert_edges_2_centers([-180:18:180]));

countsECDG=sumAngleCounts{1};
countsDGCA3=sumAngleCounts{2};
countsCA3CA1=sumAngleCounts{3};
countsECCA3=sumAngleCounts{5};
all_counts=[{countsECDG},{countsDGCA3},{countsCA3CA1},{countsECCA3}];

anglesECDG=repelem(bin_centers,countsECDG);
anglesDGCA3=repelem(bin_centers,countsDGCA3);
anglesCA3CA1=repelem(bin_centers,countsCA3CA1);
anglesECCA3=repelem(bin_centers,countsECCA3);
all_angles=[{anglesECDG},{anglesDGCA3},{anglesCA3CA1},{anglesECCA3}];

myPairs=nchoosek(1:4,2);
myPvals=myPairs;
for nPairs=1:size(myPairs)
    idx=[ones(1,sum(all_counts{myPairs(nPairs,1)})),2*ones(1,sum(all_counts{myPairs(nPairs,2)}))];

    % [pval,stats]=circ_wwtest([all_angles{myPairs(nPairs,1)},all_angles{myPairs(nPairs,2)}],idx);
    [pval,stats]=circ_ktest(all_angles{myPairs(nPairs,1)},all_angles{myPairs(nPairs,2)});

    myPvals(nPairs,3)=pval;
end