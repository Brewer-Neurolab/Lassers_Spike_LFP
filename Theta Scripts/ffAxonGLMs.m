%% GLM for Feedforward Axons

%% Setup
clear
clc

parent_axons_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

parent_wells_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes";
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
        if ~isempty(axon_spikes{fi}.up_ff)
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
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat")

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;
%% Compute GLM

clc

wellElecs=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\matching_table_wells_CW.mat");
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

    %filter for theta
    [A,B,C,D]=butter(8,10/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,4/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    theta=filtfilt(sos,g,LFP_filt);

    %filter for delta
    [A,B,C,D]=butter(8,4/(re_fs/2),'low');
    [sos,g]=ss2sos(A,B,C,D);
    LFP_filt=filtfilt(sos,g,data);

    [A,B,C,D]=butter(8,0.5/(re_fs/2),'high');
    [sos,g]=ss2sos(A,B,C,D);
    delta=filtfilt(sos,g,LFP_filt);

    %get axon spikes
    % axon_spikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
    % axon_spikes=axon_spikes.cluster_class(:,2)/1000;
    % axon_spike_train=zeros(1,length(t));
    % axon_spike_train(ismembertol(t,axon_spikes))=1;

    % GLM on theta regions deletion
    % plot axon data tagged
    thresh_mult=1;

    %define max number of samples for combining LFPs
    nsamples_combine_thresh=(1/10)*re_fs*3;
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest theta cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/10)*minLFPCycles*re_fs;

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(theta,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    validLFPIndex=[];
    for nEndPts=1:size(LFPEndPts,1)
        validLFPIndex=[validLFPIndex,LFPEndPts(nEndPts,1)-250:LFPEndPts(nEndPts,2)+250];
    end
    logicalValidLFPs=zeros(1,length(re_t));
    logicalValidLFPs(validLFPIndex(validLFPIndex>0 & validLFPIndex<length(t)))=1;
    logicalValidLFPs=logical(logicalValidLFPs);

    ThetaHilbert=hilbert(theta);
    ThetaAmp=abs(ThetaHilbert);
    ThetaAmp(~logicalValidLFPs)=[];
    ThetaAngle=angle(ThetaHilbert);
    ThetaAngle(~logicalValidLFPs)=[];

    DeltaHilbert=hilbert(delta);
    DeltaAmp=abs(DeltaHilbert);
    DeltaAmp(~logicalValidLFPs)=[];
    DeltaAngle=angle(DeltaHilbert);
    DeltaAngle(~logicalValidLFPs)=[];

    % targetElecs=well_spike_dyn.channel_name(well_spike_dyn.fi==ff_axon_tbl.fi(nFF) & well_spike_dyn.regi==ff_axon_tbl.subi(nFF));
    % sourceLFP_targetSpike_relations(t,re_t,logicalValidLFPs,LFPAmplitude,LFPAngles,ff_axon_tbl.fi(nFF),ff_axon_tbl.Electrode(nFF),targetElecs,well_spike_dyn,20,thresh_mult,...
    %     "D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\",...
    %     "C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Images\Theta")
    % 
    % close all

    myWellElecs=wellElecs.electrode(wellElecs.subregion==ff_axon_tbl.FFReg(nFF));
    
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
        glmTbl.ThetaAmp=ThetaAmp';
        % glmTbl.ThetaAngle=discretize(ThetaAngle',angle_edges,angle_edges(2:end));
        glmTbl.ThetaAngle=ThetaAngle';
        glmTbl.SinThetaAngle=sin(glmTbl.ThetaAngle);
        glmTbl.CosThetaAngle=cos(glmTbl.ThetaAngle);
        glmTbl.DeltaAmp=DeltaAmp';
        % glmTbl.DeltaAngle=discretize(resample(DeltaAngle',re_fs,fs),angle_edges,diff(angle_edges));
        glmTbl.DeltaAngle=DeltaAngle';
        glmTbl.SinDeltaAngle=sin(glmTbl.DeltaAngle);
        glmTbl.CosDeltaAngle=cos(glmTbl.DeltaAngle);

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
        modelspec='WellSpikes ~ ThetaAmp*(ThetaAngle+SinThetaAngle+CosThetaAngle)';
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

%% plot log10(p) for cos v amp
% unique_sources=unique(glmTblAll(:,[1,2,3]));

% for nAxons=1:height(unique_sources)
amp=[];
interaction=[];
colorVar=string();
for nConnect=1:height(glmTblAll)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    amp(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(2));
    interaction(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(8));
    colorVar(nConnect)=string(glmTblAll.fi(nConnect))+glmTblAll.source_elec(nConnect);
end

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),40)
end
hold off

xlabel("-log10 p Amp")
ylabel("-log10 p Amp:Cos")

ax=gca;
ax.XScale="log";
ax.YScale="log";

ylim([0.001,50])
xlim([0.001,500])

xline(-log10(0.0009))
yline(-log10(0.0009))

legend([unique_sources,'',''],'Location','southeast')
% ax.tick
% end
% hold off
nsig=sum(amp>-log10(0.0009) & interaction>-log10(0.0009));

title("#Significant Connections: "+nsig)
set(gca,"FontSize",20)
%% plot log10(p) for amp v pval
% unique_sources=unique(glmTblAll(:,[1,2,3]));

% for nAxons=1:height(unique_sources)
amp=[];
interaction=[];
colorVar=string();
for nConnect=1:height(glmTblAll)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    amp(nConnect)=-log10(glmTblAll.mdlPVal(nConnect));
    interaction(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(2));
    colorVar(nConnect)=string(glmTblAll.fi(nConnect))+glmTblAll.source_elec(nConnect);
end

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),40)
end
hold off

xlabel("-log10 p mdl")
ylabel("-log10 p Amp")

ax=gca;
ax.XScale="log";
ax.YScale="log";

ylim([0.001,50])
xlim([0.001,500])

xline(-log10(0.0009))
yline(-log10(0.0009))

legend([unique_sources,'',''],'Location','southeast')
% ax.tick
% end
% hold off
nsig=sum(amp>-log10(0.0009) & interaction>-log10(0.0009));

title("#Significant Connections: "+nsig)
set(gca,"FontSize",20)
%% plot log10(p) for amp:cos v pval EC-DG
% unique_sources=unique(glmTblAll(:,[1,2,3]));

% for nAxons=1:height(unique_sources)
amp=[];
interaction=[];
colorVar=string();
glmTblECDG=glmTblAll(glmTblAll.source_reg=="EC-DG",:);
for nConnect=1:height(glmTblECDG)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    amp(nConnect)=-log10(glmTblECDG.mdl{nConnect}.Coefficients.pValue(8));
    interaction(nConnect)=-log10(glmTblECDG.mdl{nConnect}.Coefficients.pValue(2));
    colorVar(nConnect)=string(glmTblECDG.fi(nConnect))+glmTblECDG.source_elec(nConnect);
end

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    thisColor=cVec(colorVar==unique_sources(nAxons),:);
    scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),40,"MarkerEdgeColor",thisColor(1,:),"LineWidth",1.5)
end
hold off

xlabel("-log10 p model")
ylabel("-log10 p Amp:cos Angle")

ax=gca;
ax.XScale="log";
ax.YScale="log";

% ylim([0.001,50])
% xlim([0.001,500])

BonferroniP=0.05/height(glmTblECDG);

xline(-log10(BonferroniP))
yline(-log10(BonferroniP))

legend([unique_sources,'',''],'Location','southeast')
% ax.tick
% end
% hold off
nsig=sum(amp>-log10(BonferroniP) & interaction>-log10(BonferroniP));

title("#Significant Connections: "+nsig+"/"+height(glmTblECDG))
set(gca,"FontSize",20)

axis square
%% plot log10(p) for amp:sin v pval
% unique_sources=unique(glmTblAll(:,[1,2,3]));

% for nAxons=1:height(unique_sources)
amp=[];
interaction=[];
colorVar=string();
for nConnect=1:height(glmTblAll)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    amp(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(7));
    interaction(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(2));
    colorVar(nConnect)=string(glmTblAll.fi(nConnect))+glmTblAll.source_elec(nConnect);
end

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),40)
end
hold off

xlabel("-log10 p mdl")
ylabel("-log10 p Amp:cos Angle")

ax=gca;
ax.XScale="log";
ax.YScale="log";

ylim([0.001,50])
xlim([0.001,500])

xline(-log10(0.0009))
yline(-log10(0.0009))

legend([unique_sources,'',''],'Location','southeast')
% ax.tick
% end
% hold off
nsig=sum(amp>-log10(0.0009) & interaction>-log10(0.0009));

title("#Significant Connections: "+nsig)
set(gca,"FontSize",20)
%% plot log10(p) for amp:angle v pval
% unique_sources=unique(glmTblAll(:,[1,2,3]));

% for nAxons=1:height(unique_sources)
amp=[];
interaction=[];
colorVar=string();
for nConnect=1:height(glmTblAll)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    amp(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(7));
    interaction(nConnect)=-log10(glmTblAll.mdl{nConnect}.Coefficients.pValue(2));
    colorVar(nConnect)=string(glmTblAll.fi(nConnect))+glmTblAll.source_elec(nConnect);
end

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),40)
end
hold off

xlabel("-log10 p mdl")
ylabel("-log10 p Amp:cos Angle")

ax=gca;
ax.XScale="log";
ax.YScale="log";

ylim([0.001,50])
xlim([0.001,500])

xline(-log10(0.0009))
yline(-log10(0.0009))

legend([unique_sources,'',''],'Location','southeast')
% ax.tick
% end
% hold off
nsig=sum(amp>-log10(0.0009) & interaction>-log10(0.0009));

title("#Significant Connections: "+nsig)
set(gca,"FontSize",20)