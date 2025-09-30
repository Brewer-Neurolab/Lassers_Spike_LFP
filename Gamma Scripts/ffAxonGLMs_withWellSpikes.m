%% GLM for Feedforward Axons
% Uses only the wells with spikes, not all wells
%% Setup
clear
clc
close all

parent_axons_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A";
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
well_spike_dyn=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

all_reg=[interRegions,"EC-CA1"];
%% Compute GLM

clc

wellElecs=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes 5SD Min\matching_table_wells_CW.mat");
wellElecs=wellElecs.matching_table;

glmTblAll=table();
row=1;

for nFF=1:height(ff_axon_tbl)%[95,96,117]
    % Axon Setup
    data=load(fullfile(parent_axons_dir,axons_folders(ff_axon_tbl.fi(nFF)),ff_axon_tbl.Electrode(nFF)+".mat"));
    data=data.data;
    fs=25000;
    t_rec=300;
    re_fs=1000;
    re_t=0:1/re_fs:t_rec-(1/re_fs);

    %downsample data
    data=resample(data,re_fs,fs);

    %create full sampled time steps
    t=0:1/fs:t_rec-(1/fs);

    %filter for high gamma
    [A,B,C,D]=butter(8,300/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,100/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    high_gamma=filtfilt(sos,g,LFP_filt);

    %filter for low gamma
    [A,B,C,D]=butter(8,100/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,30/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    low_gamma=filtfilt(sos,g,LFP_filt);

    %get axon spikes
    % axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
    % axon_spikes=axon_spikes.cluster_class(:,2)/1000;
    % axon_spike_train=zeros(1,length(t));
    % axon_spike_train(ismembertol(t,axon_spikes))=1;

    % GLM on theta regions deletion
    % plot axon data tagged
    thresh_mult=1;

    %define max number of samples for combining LFPs
    nsamples_combine_thresh=(1/300)*re_fs*3;
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest gamma cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/300)*minLFPCycles*re_fs;

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(high_gamma,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    % LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));
    LFPAngles=angle(LFPHilbert)*(180/pi);

    validLFPIndex=[];
    for nEndPts=1:size(LFPEndPts,1)
        validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2)];
    end
    validLFPIndex=unique(validLFPIndex);
    validLFPIndex(validLFPIndex<=0 | validLFPIndex>length(re_t))=[];
    logicalValidLFPs=zeros(1,length(re_t));
    logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(re_t)))=1;
    logicalValidLFPs=logical(logicalValidLFPs);

    HighGammaHilbert=hilbert(high_gamma);
    HighGammaAmp=abs(HighGammaHilbert);
    HighGammaAmp(~logicalValidLFPs)=[];
    HighGammaAngle=angle(HighGammaHilbert);
    HighGammaAngleWhole=HighGammaAngle;
    HighGammaAngle(~logicalValidLFPs)=[];

    % Count theta cycles
    % thetaCycleCounter=[];
    % for nEndPts=1:size(LFPEndPts,1)
    %     pkRange=LFPEndPts(nEndPts,1):LFPEndPts(nEndPts,2);
    %     pkRange(pkRange<=0 | pkRange>length(re_t))=[];
    %     [~,pkLocs]=findpeaks(ThetaAngleWhole(pkRange),"MinPeakHeight",3);
    %     thisCycle=[];
    %     for nPks=1:length(pkLocs)+1
    %         if nPks==1
    %             thisCycle=[repmat(nPks,pkLocs(nPks),1)];
    %         elseif nPks==length(pkLocs)+1
    %             thisCycle=[thisCycle;repmat(nPks,LFPEndPts(nEndPts,2)-LFPEndPts(nEndPts,1)-pkLocs(nPks-1),1)];
    %         else
    %             thisCycle=[thisCycle;repmat(nPks,pkLocs(nPks)-pkLocs(nPks-1),1)];
    %         end
    %     end
    %     thetaCycleCounter=[thetaCycleCounter;thisCycle];
    % end

    LowGammaHilbert=hilbert(low_gamma);
    LowGammaAmp=abs(LowGammaHilbert);
    LowGammaAmp(~logicalValidLFPs)=[];
    LowGammaAngle=angle(LowGammaHilbert);
    LowGammaAngle(~logicalValidLFPs)=[];

    % targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF));
    % sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,ff_axon_tbl.fi(nFF),ff_axon_tbl.Electrode(nFF),targetElecs,well_spike_dyn,20,thresh_mult,...
    %     "D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\",...
    %     "C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Images\Theta")
    % 
    % close all

    % myWellElecs=wellElecs.electrode(wellElecs.subregion==ff_axon_tbl.FFReg(nFF));
    myWellElecs=well_spike_dyn.channel_name(well_spike_dyn.regi==ff_axon_tbl.subi(nFF) & well_spike_dyn.fi==ff_axon_tbl.fi(nFF));
    
    for nWell=1:length(myWellElecs)
        angle_edges=[-pi:pi/180*18:pi];
        % spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\A8_spikes.mat");
        spikes=load(fullfile(parent_wells_dir,wells_folders(ff_axon_tbl.fi(nFF)),myWellElecs(nWell)+"_spikes.mat"));
        spikes=spikes.index/1000;
        spike_train=zeros(1,length(t));
        spike_train(ismembertol(t,spikes))=1;
        glmTbl=table();
        spike_idx=round(remap(find(spike_train),1,length(t),1,300*re_fs));
        spike_train=zeros(1,300*re_fs);
        spike_train(spike_idx)=1;
        spike_train(~logicalValidLFPs)=[];

        glmTbl.WellSpikes=logical(spike_train');
        % glmTbl.ThetaAmp=ThetaAmp';
        glmTbl.HighGammaAmp=zscore(HighGammaAmp'); %toggle to zscore
        % glmTbl.ThetaAngle=discretize(ThetaAngle',angle_edges,angle_edges(2:end));
        % glmTbl.ThetaAngle=ThetaAngle';
        % glmTbl.SinThetaAngle=sin(glmTbl.ThetaAngle);
        % glmTbl.CosThetaAngle=cos(glmTbl.ThetaAngle);

        glmTbl.HighGammaAngle=zscore(HighGammaAngle'); %toggle to zscore
        glmTbl.SinHighGammaAngle=zscore(sin(glmTbl.HighGammaAngle)); %toggle to zscore
        glmTbl.CosHighGammaAngle=zscore(cos(glmTbl.HighGammaAngle)); %toggle to zscore
        
        % glmTbl.DeltaAmp=DeltaAmp';
        % glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
        % glmTbl.DeltaAngle=DeltaAngle';
        % glmTbl.SinDeltaAngle=sin(glmTbl.DeltaAngle);
        % glmTbl.CosDeltaAngle=cos(glmTbl.DeltaAngle);

        % modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
        %     'SinThetaAngle + CosThetaAngle + SinDeltaAngle + CosDeltaAngle+'...
        %     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle*'...
        %     'SinThetaAngle*CosThetaAngle*SinDeltaAngle*CosDeltaAngle'];
        % modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
        %     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
        % modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + '...
        %     'ThetaAmp*ThetaAngle'];
        % modelspec=['WellSpikes ~ ThetaAmp + ThetaAngle + DeltaAmp + DeltaAngle + '...
        %     'ThetaAmp*ThetaAngle*DeltaAmp*DeltaAngle'];
        % modelspec='WellSpikes ~ ThetaAmp*(ThetaAngle+SinThetaAngle+CosThetaAngle)';
        modelspec='WellSpikes ~ HighGammaAmp*(CosHighGammaAngle)'; %refined from all three angle variations
        % modelspec='WellSpikes ~ ThetaAmp:ThetaAngle';
        % modelspec=['WellSpikes ~ ThetaAmp*ThetaAngle*SinThetaAngle*CosThetaAngle' ...
        %     '*DeltaAmp*DeltaAngle*SinDeltaAngle*CosDeltaAngle'];

        mdl_del=fitglm(glmTbl,modelspec,'Distribution','binomial','Link','logit','LikelihoodPenalty','jeffreys-prior');

        glmTblAll.fi(row)=ff_axon_tbl.fi(nFF);
        glmTblAll.source_reg(row)=ff_axon_tbl.Subregion(nFF);
        glmTblAll.source_elec(row)=ff_axon_tbl.Electrode(nFF);
        glmTblAll.target_reg(row)=ff_axon_tbl.FFReg(nFF);
        glmTblAll.target_elec(row)=myWellElecs(nWell);
        glmTblAll.mdl{row}=mdl_del;
        % glmTblAll.mdlPVal(row)=mdl_del.devianceTest.pValue(2);
        glmTblAll.mdlPVal(row)=coefTest(mdl_del);


        row=row+1;

        disp(ff_axon_tbl.fi(nFF)+" "+myWellElecs(nWell))
    end

end

%% Create low memory glm table

for nConnect=1:height(glmTblAll)
    lowInfoGLM=struct(glmTblAll.mdl{nConnect});
    myFields=fieldnames(lowInfoGLM);
    toRemove=myFields(ismember(myFields,["Design","design_r","y_r","w_r","ObservationInfo","Data","Variables","Offset","IRLSWeights","Residuals","Leverage"]));
    glmTblAll.mdl{nConnect}=rmfield(lowInfoGLM,toRemove);
end

%% Count significant by subregion
for subi=1:length(all_reg)
    subTbl=glmTblAll(glmTblAll.source_reg==all_reg(subi),:);
    BonferroniP=0.05/height(subTbl);
    disp("#significant "+all_reg(subi)+"="+sum(subTbl.mdlPVal<BonferroniP)+"/"+height(subTbl))
end
%% Amp V Mdl by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,"mdl",2,all_reg(subi),0.05)
end

%% CosAngle V Mdl by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,"mdl",3,all_reg(subi),0.05)
end
%% CosAngle V SinAngle by subregion

% for subi=1:length(all_reg)
%     glmScatter(glmTblAll,4,5,all_reg(subi),0.05)
% end
%% Angle V Mdl by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,"mdl",3,all_reg(subi),0.05)
end
%% Cos Angle V Amp by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,2,3,all_reg(subi),0.05)
end
%% Angle V Amp by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,2,3,all_reg(subi),0.05)
end
%% Amp:CosAngle V Mld by subregion

for subi=1:length(all_reg)
    glmScatter(glmTblAll,"mdl",4,all_reg(subi),0.05)
end