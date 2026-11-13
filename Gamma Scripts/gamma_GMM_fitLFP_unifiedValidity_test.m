%% Setup
clear
clc
close all

saveDir="C:\BrewerLabResearch\GitHub\Lassers_Spike_LFP\Gamma Scripts";

parent_axons_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\downsampled tunnels\ASR\Low_Gamma"; %Change b/w low and high gamma
axons_dir=dir(parent_axons_dir);
axons_folders=string({axons_dir.name});
axons_folders=axons_folders([axons_dir.isdir]);
axons_folders=axons_folders(3:end);

figOutDir = "C:\BrewerLabResearch\Presentations&Figures\GaussFitFigsASR_LowGamma"; %Change b/w low and high gamma
if ~isfolder(figOutDir), mkdir(figOutDir); end
% ==========================================================================

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
        % if ~isempty(axon_spikes{fi}.up_ff{nelec}) | (~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC")
        if ~isempty(axon_spikes{fi}.up_ff{nelec}) | ~isempty(axon_spikes{fi}.up_fb{nelec})
            ff_axon_tbl.fi(row)=fi;
            ff_axon_tbl.Subregion(row)=axon_spikes{fi}.Subregion(nelec);
            ff_axon_tbl.interRegi(row)=find(interRegions==axon_spikes{fi}.Subregion(nelec));
            electrodes=axon_spikes{fi}.("Electrode Pairs")(nelec);
            electrodes=split(electrodes,{'-'});
            ff_axon_tbl.Electrode(row)=electrodes(1);
            tunnelReg=split(axon_spikes{fi}.Subregion(nelec),{'-'});
            ff_axon_tbl.FFReg(row)=tunnelReg(2);
            ff_axon_tbl.subi(row)=find(subregions==tunnelReg(2));

            % if isempty(axon_spikes{fi}.up_ff{nelec}) & ~isempty(axon_spikes{fi}.up_fb{nelec}) & axon_spikes{fi}.Subregion(nelec)=="CA1-EC"
            %     ff_axon_tbl.Subregion(row)="EC-CA1";
            %     ff_axon_tbl.FFReg(row)="CA1";
            %     ff_axon_tbl.subi(row)=4;
            % end

            row=row+1;
        end
    end
end

%3.5 SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

% all_reg=[interRegions,"EC-CA1"];
all_reg=[interRegions];
%all_reg=["EC-DG"];

%% Get hilbert amps of each axon

clc

%testing range
testing_idx=[find(ff_axon_tbl.Subregion=="CA3-CA1")]';

LFPTable=ff_axon_tbl;

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
    nsamples_combine_thresh=(1/100)*re_fs*3; %Change b/w low and high gamma
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest gamma cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/100)*minLFPCycles*re_fs; %Change b/w low and high gamma

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    % ===== EDIT 2: strip zero / non-finite envelope samples. log10(0) = -Inf
    %               breaks logspace, histcounts and fitgmdist. Harmless for
    %               theta, but gamma envelopes do hit exact zeros.
    LFPAmplitude = LFPAmplitude(isfinite(LFPAmplitude) & LFPAmplitude > 0);
    % Keep every LFP cell in the same orientation.
    LFPAmplitude = LFPAmplitude(:);
    % ======================================================================

    LFPTable.LFP{nFF}=LFPAmplitude;
    disp(nFF+"/"+height(ff_axon_tbl))
end

% Identify rows containing usable LFP amplitude data.
hasLFP = false(height(LFPTable),1);

for r = 1:height(LFPTable)
    hasLFP(r) = ...
        ~isempty(LFPTable.LFP{r}) && ...
        isnumeric(LFPTable.LFP{r}) && ...
        any(isfinite(LFPTable.LFP{r}) & LFPTable.LFP{r} > 0);
end

processedRows = find(hasLFP);

fprintf('Usable LFP data: %d/%d electrode pairs\n', ...
    numel(processedRows),height(LFPTable));

if isempty(processedRows)
    error('No usable LFP amplitude data were extracted.');
end

% Create one robust global amplitude range for later sections.
allAmplitudes = vertcat(LFPTable.LFP{processedRows});
allAmplitudes = allAmplitudes( ...
    isfinite(allAmplitudes) & allAmplitudes > 0);

globalAmpRange = [min(allAmplitudes),max(allAmplitudes)];

%% get average histograms of each and fit gaussian
close all
avgAmpOverThresh=[];
sdAmpOverThresh=[];
lfpPropsTbl=table();
row=1;
axonsCounted=0;
fitVec=[];
for i=1:length(all_reg)
    figure
    regRows = find(LFPTable.Subregion == all_reg(i) & hasLFP);
    if isempty(regRows)
        warning('%s has no usable LFP data; skipping.',all_reg(i));
        continue
    end
    axonsCounted=axonsCounted+length(regRows);
    disp(axonsCounted+"/"+height(LFPTable))
    %ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    regionAmplitudes = vertcat(LFPTable.LFP{regRows});
    regionAmplitudes = regionAmplitudes( ...
        isfinite(regionAmplitudes) & regionAmplitudes > 0);

    ampRange = [min(regionAmplitudes),max(regionAmplitudes)];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    ampPts=[];
    cdfPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        % OLD: figure; pts=histogram(LFPTable.LFP{regRows(naxons)},binEdges);
        % figure; mycdf=histogram(...,"Normalization","cdf");
        ampPts(naxons,:) = histcounts(LFPTable.LFP{regRows(naxons)}, binEdges);
        cdfPts(naxons,:) = histcounts(LFPTable.LFP{regRows(naxons)}, binEdges, 'Normalization','cdf');
        % ampPts(naxons,:)=pts.BinCounts;
        % cdfPts(naxons,:)=mycdf.Values;
        axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));
    end

    myColors=distinguishable_colors(size(ampPts,1),{'w','r','k'});
    % myColors=[myColors,repmat(0.5,[size(myColors,1),1])];

    plot(10.^convert_edges_2_centers(log10(binEdges)),ampPts,'LineWidth',1)
    ax1=gca;
    ax1.XScale="log";
    ax1.ColorOrder=myColors;
    xlabel("Amplitude (uV)")
    ylabel("Count")

    meanPts=mean(ampPts,1);
    sdPts=std(ampPts,[],1);
    % plot(convert_edges_2_centers(log(binEdges)),meanPts,'k')
    % hold on
    % errorbar(convert_edges_2_centers(log(binEdges)),meanPts,sdPts)
    f=fit(convert_edges_2_centers(log(binEdges))',meanPts','gauss1');
    fitvals=f(convert_edges_2_centers(log(binEdges))');
    fitVec{i}=f;

    % ===== EDIT 3: stash the per-region bin edges + fit values so the new
    %               electrode-axon-pair figure cell can reproduce these curves
    %               without recomputing the fit.
    regBinEdges{i} = binEdges;
    regFitVals{i}  = fitvals;
    % ======================================================================

    figure
    plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals,'r--')
    hold on

    plot(exp(convert_edges_2_centers(log(binEdges)))',meanPts,'k-.')
    errorbar(exp(convert_edges_2_centers(log(binEdges)))',meanPts,sdPts,'Color','k','LineStyle','none')
    hold off

    ax2=gca;
    ax2.XScale="log";
    xlim(ampRange)

    mu=f.b1;
    sigma=f.c1/sqrt(2);

    m=exp(mu+(sigma^2/2));
    v=exp(2*mu+sigma^2)*(exp(sigma^2)-1);
    ampThreshold = m+2*sqrt(v);
    xline(ax1,ampThreshold)
    legend(ax1,axonNames)
    xline(ax2,ampThreshold)

    % ===== EDIT 4: clamp last_bin. If ampThreshold lands past the final edge (common
    %               in high gamma, short tail) last_bin(1) errored on empty.
    % last_bin=find(binEdges>ampThreshold);
    % last_bin=last_bin(1);
    last_bin=find(binEdges>ampThreshold);
    if isempty(last_bin)
        last_bin = numel(binEdges)-1;
    else
        last_bin = last_bin(1);
    end
    last_bin = min(last_bin, size(ampPts,2));
    % ======================================================================

    title(ax1,all_reg(i)+" "+sum(any(ampPts(:,last_bin:end),2))+"/"+size(ampPts,1))

    figure
    hold on
    % plot(10.^convert_edges_2_centers(log10(binEdges)),ampPts,'LineWidth',4,'Color',[0,0,0,0.5])
    for ncolors=1:size(myColors,1)
        plot(10.^convert_edges_2_centers(log10(binEdges)),ampPts(ncolors,:),'LineWidth',4,'Color',[myColors(ncolors,:)])
    end

    plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals,'r--','LineWidth',6)
    % plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals*1.02,'k','LineWidth',2)
    % plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals*0.98,'k','LineWidth',2)
    area(exp(convert_edges_2_centers(log(binEdges)))',fitvals,'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none')
    xlabel("Amplitude (uV)")
    ylabel("Count")
    hold off

    ax3=gca;
    ax3.XScale="log";
    % ax3.ColorOrder=myColors;
    xlim([-Inf,1.5*10^3])
    xl=xlim;
    yl=ylim;
    xticks(logspace(-3,4,8))
    xticklabels("10^{"+[-3:4]+"}")
    % ylim([0,10^3])
    axis square
    l=legend(ax3,[axonNames,"Gauss Fit"]);
    l.AutoUpdate="off";
    xline(ax3,ampThreshold,'LineWidth',4,'Color','k','LineStyle','--')
    % title(ax3,all_reg(i)+" "+sum(any(ampPts(:,last_bin:end),2))+"/"+size(ampPts,1))
    % title(ax3,all_reg(i)+" "+size(ampPts,1))
    t=text(0.025,0.95,ax3,all_reg(i)+" "+"n="+size(ampPts,1),'Units','normalized','FontSize',30);
    set(ax3,"FontSize",28)
    box on
    ax3.LineWidth=4;
    yticks(linspace(0,30000,16))

    ampWeightCount=ampPts.*10.^convert_edges_2_centers(log10(binEdges));
    %
    % avgAmpOverThresh(i)=(mean(ampPts(any(ampPts(:,last_bin:end),2),last_bin:end),"all"));
    % sdAmpOverThresh=(std(ampPts(any(ampPts(:,last_bin:end),2),last_bin:end),[],"all"));
    axonsAboveThresh=find(any(ampPts(:,last_bin:end),2));
    for naxons=1:length(axonsAboveThresh)
        lfpPropsTbl.Subregion(row)=all_reg(i);
        lfpPropsTbl.AxonName(row)=axonNames(axonsAboveThresh(naxons));
        lfpPropsTbl.countIntegration(row)=sum(ampPts(axonsAboveThresh(naxons),last_bin:end));
        lfpPropsTbl.sIntegration(row)=sum(ampPts(axonsAboveThresh(naxons),last_bin:end))/re_fs;
        lfpPropsTbl.ampIntegration(row)=sum(ampWeightCount(axonsAboveThresh(naxons),last_bin:end));
        row=row+1;
    end

    figure
    hold on
    % plot(10.^convert_edges_2_centers(log10(binEdges)),ampPts,'LineWidth',4,'Color',[0,0,0,0.5])
    for ncolors=1:size(myColors,1)
        plot(10.^convert_edges_2_centers(log10(binEdges)),cdfPts(ncolors,:),'LineWidth',4,'Color',[myColors(ncolors,:)])
    end

    % regioncdf=histcounts([LFPTable.LFP{regRows}],binEdges,"Normalization","cdf");
    regionAmplitudes = vertcat(LFPTable.LFP{regRows});
    regioncdf = histcounts(regionAmplitudes,binEdges, "Normalization","cdf");
    plot(10.^convert_edges_2_centers(log10(binEdges)),regioncdf,'r--','LineWidth',8)

    % plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals,'r--','LineWidth',6)
    % plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals*1.02,'k','LineWidth',2)
    % plot(exp(convert_edges_2_centers(log(binEdges)))',fitvals*0.98,'k','LineWidth',2)
    % area(exp(convert_edges_2_centers(log(binEdges)))',fitvals,'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none')
    xlabel("Amplitude (uV)")
    ylabel("Count")
    hold off

    ax4=gca;
    ax4.XScale="log";
    % ax4.YScale="log";
    % ax3.ColorOrder=myColors;
    xlim([-Inf,1.5*10^3])
    ylim([0.1,Inf])
    xl=xlim;
    yl=ylim;
    xticks(logspace(-3,4,8))
    xticklabels("10^{"+[-3:4]+"}")
    % ylim([0,10^3])
    axis square
    l=legend(ax4,axonNames);
    l.AutoUpdate="off";
    xline(ax4,ampThreshold,'LineWidth',4,'Color','k','LineStyle','--')
    % title(ax3,all_reg(i)+" "+sum(any(ampPts(:,last_bin:end),2))+"/"+size(ampPts,1))
    % title(ax3,all_reg(i)+" "+size(ampPts,1))
    t=text(0.025,0.95,ax4,all_reg(i)+" "+"n="+size(ampPts,1),'Units','normalized','FontSize',30);
    set(ax4,"FontSize",28)
    box on
    ax4.LineWidth=4;
    yticks(linspace(0,30000,16))

end

% calculate means
ampMeanTab=table();
for i=1:length(all_reg)
    ampMeanTab.subregion(i)=all_reg(i);
    ampMeanTab.mean(i)=mean(log10(lfpPropsTbl.ampIntegration(lfpPropsTbl.Subregion==all_reg(i))));
    ampMeanTab.sd(i)=std(log10(lfpPropsTbl.ampIntegration(lfpPropsTbl.Subregion==all_reg(i))));
    ampMeanTab.min(i)=ampMeanTab.mean(i)-(2*ampMeanTab.sd(i));
    ampMeanTab.max(i)=ampMeanTab.mean(i)+(2*ampMeanTab.sd(i));
    ampMeanTab.n(i)=sum(lfpPropsTbl.Subregion==all_reg(i));
    ampMeanTab.meanCount(i)=mean(lfpPropsTbl.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
    %remove outliers +/- 2SD
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & (lfpPropsTbl.ampIntegration<=10^ampMeanTab.min(i) | lfpPropsTbl.ampIntegration>=10^ampMeanTab.max(i)),:)=[];
    % remove integrations less than 50000
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.ampIntegration<=50000,:)=[];
    % remove time less than 5 seconds
    lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.sIntegration<=30,:)=[];
end

%new means
newAmpMeanTab=table();
for i=1:length(all_reg)
    newAmpMeanTab.subregion(i)=all_reg(i);
    newAmpMeanTab.mean(i)=mean(log10(lfpPropsTbl.ampIntegration(lfpPropsTbl.Subregion==all_reg(i))));
    newAmpMeanTab.sd(i)=std(log10(lfpPropsTbl.ampIntegration(lfpPropsTbl.Subregion==all_reg(i))));
    newAmpMeanTab.min(i)=newAmpMeanTab.mean(i)-(2*newAmpMeanTab.sd(i));
    newAmpMeanTab.max(i)=newAmpMeanTab.mean(i)+(2*newAmpMeanTab.sd(i));
    newAmpMeanTab.n(i)=sum(lfpPropsTbl.Subregion==all_reg(i));
    newAmpMeanTab.meanCount(i)=mean(lfpPropsTbl.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
end

%
regCats=categorical(all_reg);
regCats=reordercats(regCats,all_reg);
figure
bar(regCats,newAmpMeanTab.mean)
hold on
errorbar(regCats,newAmpMeanTab.mean,newAmpMeanTab.sd./sqrt(newAmpMeanTab.n),'LineStyle','none','Color','k')
hold off
% ===== EDIT 5: hard-coded ylim was tuned for theta and hides gamma data.
% ylim([5,7])
ylim('auto')
% ======================================================================
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)

% %% Stats calc
% y=log10(lfpPropsTbl.ampIntegration);
% g={lfpPropsTbl.Subregion};
% 
% [p,tbl,stats]=anovan(y,g);
% [c,m]=multcompare(stats,'CriticalValueType','hsd');
% 
% figure
% bar(regCats,m(:,1),'FaceColor','#75ebeb')
% hold on
% errorbar(regCats,m(:,1),m(:,2),'LineStyle','none','Color','k','LineWidth',3)
% hold off
% ylim([5,8])
% % ylabel("Log_{10} Amplitude Integration")
% ylabel("Mean Amplitude Integral x10^6")
% yticks(0:1:10)
% yticklabels(10.^(0:1:10)./10^6)
% axis square
% set(gca,"FontSize",28)
% ax=gca;
% ax.YAxis.MinorTick="on";
% ax.YAxis.MinorTickValues=0:0.2:10;
% box on
% ax.LineWidth=4;
% ax.TickLength=[0.05,0.025];

%% Gaussian Mixture Model  (REPLACES the old GMM section)
close all
rng('default')

% ---------- Settings: shared by fitting, figures AND final table ----------
P = struct();
P.minMu2_uV       = 30;     % absolute floor for mu2 -- band dependent!
P.minMu2_log      = log10(P.minMu2_uV);
P.minWeight2      = 0.01;
P.minCount2       = 50;
P.minSeparation   = 2;      % mu2 > mu1 + minSeparation*SD1
P.minDeltaBIC     = 10;     % BIC1 - BIC2 must exceed this
P.requireBoundary = true;   % require posterior crossing between mu1 and mu2
P.nReplicates     = 10;     % k-means++ restarts for the automatic K=2 fit
P.showRejectedMixture = false; % true: draw rejected K=2 mixture as thin grey line
P.clearOldFigures = false;  % true: delete old "GMM *.png" before saving
options = statset('MaxIter',1000,'TolFun',1e-8);

if P.clearOldFigures
    delete(fullfile(figOutDir,"GMM *.png"));
end

ampRange = globalAmpRange;
binEdges = logspace(log10(ampRange(1)),log10(ampRange(2)),101);
edgesLog = log10(binEdges);

% ---------- Preallocate: every row gets every column ----------
gofTable = LFPTable(:,1:6);
nRows = height(gofTable);

numVars = ["nSamples","maxAmplitude_uV","muK1","sdK1","BIC1", ...
    "NegativeLogLikelihood","BIC2","deltaBIC", ...
    "mu1","mu2","sigma1","sigma2","std1","std2","weights1","weights2", ...
    "ampUnweighted1","ampUnweighted2","ampWeighted1","ampWeighted2", ...
    "boundary","sepInSD1","ashmanD","n2Est", ...
    "muNoise","sdNoise","wNoise", ...
    "integration_log10","integration","percentIntegrated","nMu1","nMu2", ...
    "ttest_h","ttest_p","ttest_negLog10P"];
for vn = numVars, gofTable.(char(vn)) = nan(nRows,1); end

logVars = ["thresholdPass","hasBoundary","bicPass","wPass","sepPass", ...
    "mu2Valid","ttestPerformed"];
for vn = logVars, gofTable.(char(vn)) = false(nRows,1); end

gofTable.fitStatus      = repmat("not_processed",nRows,1);
gofTable.fitSource      = repmat("",nRows,1);
gofTable.failedCriteria = repmat("",nRows,1);
gofTable.minMu2_log     = repmat(P.minMu2_log,nRows,1);
gofTable.gmm    = cell(nRows,1);
gofTable.myData = cell(nRows,1);

for i = 1:numel(processedRows)
    r = processedRows(i);
    name = gofTable.Subregion(r)+" FID "+gofTable.fi(r)+" "+gofTable.Electrode(r);

    myData = log10(LFPTable.LFP{r}(:));
    myData = myData(isfinite(myData));
    N = numel(myData);
    gofTable.myData{r} = myData;
    gofTable.nSamples(r) = N;
    gofTable.maxAmplitude_uV(r) = 10^max(myData);

    if N < 2*P.minCount2
        gofTable.fitStatus(r) = "too_few_samples";
        gofTable.failedCriteria(r) = "too few samples";
        continue
    end

    regValue = max(1e-10,1e-6*var(myData,1));

    % ---- K = 1: always fitted. This IS the noise model when mu2 is rejected.
    gmm1 = fitgmdist(myData,1,'Options',options,'RegularizationValue',regValue);
    muK1 = gmm1.mu;
    sdK1 = sqrt(gmm1.Sigma);
    gofTable.muK1(r) = muK1;
    gofTable.sdK1(r) = sdK1;
    gofTable.BIC1(r) = gmm1.BIC;

    % ---- K = 2: targeted start (if possible) + automatic restarts; best NLL wins
    [gmm,src] = fitBestTwoComponent(myData,P,options,regValue);

    mu = [NaN;NaN]; vr = [NaN;NaN]; w = [NaN;NaN]; sd = [NaN;NaN]; bnd = NaN;

    if isempty(gmm)
        gofTable.fitStatus(r) = "k2_fit_failed";
        C = evaluateMu2Criteria(mu,sd,w,bnd,gmm1.BIC,NaN,N,P);
        C.failedText = "K=2 fit failed";
    else
        gofTable.fitStatus(r) = "fitted";
        gofTable.fitSource(r) = src;
        gofTable.gmm{r} = gmm;

        % Sort: component 1 = low mean, component 2 = high mean
        [mu,idx] = sort(gmm.mu(:),'ascend');
        vAll = squeeze(gmm.Sigma); vAll = vAll(:);
        wAll = gmm.ComponentProportion(:);
        vr = vAll(idx); w = wAll(idx); sd = sqrt(vr);

        bnd = findBetweenMeanBoundary(mu,vr,w,min(myData),max(myData));
        C = evaluateMu2Criteria(mu,sd,w,bnd,gmm1.BIC,gmm.BIC,N,P);

        % ---- Store ALL statistics regardless of pass/fail ----
        gofTable.NegativeLogLikelihood(r) = gmm.NegativeLogLikelihood;
        gofTable.BIC2(r)     = gmm.BIC;
        gofTable.deltaBIC(r) = gmm1.BIC - gmm.BIC;
        gofTable.mu1(r) = mu(1);     gofTable.mu2(r) = mu(2);
        gofTable.sigma1(r) = vr(1);  gofTable.sigma2(r) = vr(2);   % variances
        gofTable.std1(r) = sd(1);    gofTable.std2(r) = sd(2);
        gofTable.weights1(r) = w(1); gofTable.weights2(r) = w(2);
        ampU = 1./sqrt(2*pi*vr);
        gofTable.ampUnweighted1(r) = ampU(1);
        gofTable.ampUnweighted2(r) = ampU(2);
        gofTable.ampWeighted1(r)   = w(1)*ampU(1);
        gofTable.ampWeighted2(r)   = w(2)*ampU(2);
        gofTable.boundary(r) = bnd;
        gofTable.sepInSD1(r) = (mu(2)-mu(1))/sd(1);
        gofTable.ashmanD(r)  = sqrt(2)*abs(mu(2)-mu(1))/sqrt(vr(1)+vr(2));
        gofTable.n2Est(r)    = w(2)*N;
        gofTable.nMu1(r) = sum(abs(myData-mu(1)) < 2*sd(1));
        gofTable.nMu2(r) = sum(abs(myData-mu(2)) < 2*sd(2));

        if isfinite(bnd)
            above = myData > bnd;
            gofTable.integration_log10(r) = sum(myData(above));
            gofTable.integration(r)       = sum(10.^myData(above));
            gofTable.percentIntegrated(r) = 100*mean(above);

            % Descriptive only (groups are defined by the fitted boundary)
            if sum(above) > 1 && sum(~above) > 1
                [h,p] = ttest2(myData(~above),myData(above),'Vartype','unequal');
                gofTable.ttest_h(r) = h;
                gofTable.ttest_p(r) = p;
                gofTable.ttestPerformed(r) = true;
                if p > 0
                    gofTable.ttest_negLog10P(r) = -log10(p);
                else
                    gofTable.ttest_negLog10P(r) = Inf;
                end
            end
        end
    end

    % ---- Single source of truth for pass/fail ----
    gofTable.thresholdPass(r) = C.threshold;
    gofTable.hasBoundary(r)   = C.boundary;
    gofTable.bicPass(r)       = C.bic;
    gofTable.wPass(r)         = C.weight;
    gofTable.sepPass(r)       = C.separation;
    gofTable.mu2Valid(r)      = C.valid;
    gofTable.failedCriteria(r)= C.failedText;

    % Noise estimate used downstream: K=2 mu1 if bimodal, otherwise the K=1 fit
    if C.valid
        gofTable.muNoise(r) = mu(1); gofTable.sdNoise(r) = sd(1); gofTable.wNoise(r) = w(1);
    else
        gofTable.muNoise(r) = muK1;  gofTable.sdNoise(r) = sdK1;  gofTable.wNoise(r) = 1;
    end

    fig = plotGmmDiagnostic(myData,edgesLog,name,mu,sd,w,bnd,C,muK1,sdK1,P);
    saveas(fig,fullfile(figOutDir,"GMM muSort "+name+".png"));
    close(fig);

    if C.valid, status = "PASS"; else, status = "FAIL (" + C.failedText + ")"; end
    fprintf('%d/%d  %s : %s\n',i,numel(processedRows),char(name),char(status));
end

% ---------- Export (all rows, all statistics) ----------
xlsxTable = removevars(gofTable,{'gmm','myData'});
xlsxTable.ttest_p_display = repmat("not tested",height(xlsxTable),1);
for k = 1:height(xlsxTable)
    if xlsxTable.ttestPerformed(k)
        if xlsxTable.ttest_p(k) == 0
            xlsxTable.ttest_p_display(k) = "< numerical precision";
        else
            xlsxTable.ttest_p_display(k) = string(sprintf('%.6e',xlsxTable.ttest_p(k)));
        end
    end
end
writetable(xlsxTable,fullfile(figOutDir,"byAxonGMMFits_muSort_final.xlsx"));

%% Single axon Dunnett testing and removal for Mu2
% OLD ANALYSIS DISABLED:
% The Dunnett procedure gives each individual mu2 group only one
% observation and destructively removes fitted boundaries.
%{
if ~exist("gofTable","var")
    gofTable = readtable( ...
    fullfile(figOutDir,"byAxonGMMFits_muSort.xlsx")); 
end

%ampRange=[min(cell2mat(LFPTable.LFP),[],"all"),max(cell2mat(LFPTable.LFP),[],"all")];
ampRange = globalAmpRange;
binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
binCenters=convert_edges_2_centers(log10(binEdges))';
x_grid = linspace(log10(ampRange(1)),log10(ampRange(2)),100)';
gauss = @(x, m, s) (1 / sqrt(2*pi*s)) .* exp(-((x - m).^2) / (2*s));

% pvalVec

for i=1:length(all_reg)

    % Anova variables
    x=[];
    g=[];

    figure
    hold on
    regRows = find( ...
        LFPTable.Subregion == all_reg(i) & hasLFP);

    if isempty(regRows)
        warning('%s has no usable LFP data; skipping.',all_reg(i));
        continue
    end
    axonNames=string();
    noiseFitMat=[];
    goodMu2=[];
    weight2=[];
    axonID=[];
    axonNum=[];
    % [maxAmp,maxAmpIdx]=cellfun(@max,LFPTable.LFP(regRows));
    % [~,maxAmpAxonIdx]=max(maxAmp);
    % Obtain one scalar maximum per electrode.
    maxAmp = nan(numel(regRows),1);

    for k = 1:numel(regRows)
        thisLFP = LFPTable.LFP{regRows(k)};
        thisLFP = thisLFP(isfinite(thisLFP));

        if ~isempty(thisLFP)
            maxAmp(k) = max(thisLFP(:));
        end
    end

    if any(isfinite(maxAmp))
        [~,maxAmpAxonIdx] = max(maxAmp,[],'omitnan');
    else
        maxAmpAxonIdx = [];
    end
    for naxons=1:length(regRows)
        if iscell(gofTable.boundaries(regRows(naxons)))
            if ~isempty(gofTable.boundaries{regRows(naxons)})
                mu2=gofTable.mu2(regRows(naxons));
                myWeight=gofTable.weights2(regRows(naxons));
                sigma2=gofTable.sigma2(regRows(naxons));
                if mu2>gofTable.boundaries{regRows(naxons)}
                    goodMu2=[goodMu2,mu2];
                    weight2=[weight2,myWeight];
                    axonID=[axonID,gofTable.Subregion(regRows(naxons))+" FID "+gofTable.fi(regRows(naxons))+" "+gofTable.Electrode(regRows(naxons))+" #"+regRows(naxons)];
                    axonNum=[axonNum;regRows(naxons)];
                end
            end
        else
            if ~isempty(gofTable.boundaries(regRows(naxons)))
                mu2=gofTable.mu2(regRows(naxons));
                myWeight=gofTable.weights2(regRows(naxons));
                sigma2=gofTable.sigma2(regRows(naxons));
                if mu2>gofTable.boundaries(regRows(naxons))
                    goodMu2=[goodMu2,mu2];
                    weight2=[weight2,myWeight];
                end
            end
        end

        %if naxons==maxAmpAxonIdx
        if ~isempty(maxAmpAxonIdx) && naxons == maxAmpAxonIdx
            histogram(log10(LFPTable.LFP{regRows(naxons)}), log10(binEdges), 'Normalization', 'pdf', ...
        'FaceColor', [0.2 0.6 1.0], 'FaceAlpha', 0.5)
        end
    end
    if isempty(goodMu2)
        warning('%s: no valid mu2, skipping Dunnett.',all_reg(i));
        continue
    end
    % %myMu1=gofTable.mu1(regRows);
    % validMu1Rows = regRows( ...
    %     gofTable.fitStatus(regRows) == "fitted" & ...
    %     isfinite(gofTable.mu1(regRows)));
    % % myMu1 = gofTable.mu1(validMu1Rows);
    % % mu1Avg = mean(myMu1);
    % % newSD = std(gofTable.mu1(regRows));
    % % meanPrior1=max(gofTable.weights1(regRows));
    % validMu1Rows = regRows( ...
    %     gofTable.fitStatus(regRows) == "fitted" & ...
    %     isfinite(gofTable.mu1(regRows)));
    % myMu1 = gofTable.mu1(validMu1Rows);
    % newSD = std(myMu1);
    % meanPrior1 = max(gofTable.weights1(validMu1Rows),[],'omitnan');
    % plot(x_grid,gauss(x_grid,mu1Avg,newSD).*meanPrior1,"LineWidth",4,"Color","b")
    
    % Use only successfully fitted, finite mu1 values.
    validMu1Rows = regRows( ...
        gofTable.fitStatus(regRows) == "fitted" & ...
        isfinite(gofTable.mu1(regRows)));

    myMu1 = gofTable.mu1(validMu1Rows);

    if isempty(myMu1)
        warning('%s has no successfully fitted mu1 values; skipping Dunnett.', ...
            all_reg(i));
        close(gcf)
        continue
    end

    mu1Avg = mean(myMu1,'omitnan');
    newVariance = var(myMu1,'omitnan');
    meanPrior1 = max(gofTable.weights1(validMu1Rows),[],'omitnan');

    if isfinite(newVariance) && newVariance > 0 && isfinite(meanPrior1)
        plot(x_grid, ...
            gauss(x_grid,mu1Avg,newVariance).*meanPrior1, ...
            "LineWidth",4,"Color","b")
    else
        warning('%s: insufficient variation to plot the regional mu1 curve.', ...
            all_reg(i));
    end

    mu2Avg = mean(goodMu2,'omitnan');
    mu2Variance = var(goodMu2,'omitnan');
    meanPrior2 = max(weight2,[],'omitnan');

    % A single mu2 value has zero/undefined regional variance.
    if ~isfinite(mu2Variance) || mu2Variance <= 0
        mu2Variance = eps;
    end

    mu2Curve = ...
        gauss(x_grid,mu2Avg,mu2Variance).*meanPrior2;

    plot(x_grid,mu2Curve, ...
        "LineWidth",4,"Color","r")

    area(x_grid,mu2Curve, ...
        'FaceColor','r', ...
        'FaceAlpha',0.25, ...
        'EdgeColor','none')
    hold off

    ax=gca;
    xticks([-2:4])
    title(all_reg(i))
    xlabel('Log_{10} \muV');
    ylabel('Probability Density');
    axis square
    box on
    set(ax,"FontSize",32)
    ax.LineWidth=4;
    ax.TickLength=[0.05,0.025];
    xlim(log10(ampRange))

    ctrlName = all_reg(i)+" mu1";
    x=[x;myMu1;goodMu2'];
    g=[g;repmat(ctrlName,[length(myMu1),1]);axonID'+" mu2"];
    g = categorical(g);
    g = reordercats(g,[string(ctrlName); setdiff(string(categories(g)),string(ctrlName))]);
    figure
    [~,~,stats]=anovan(x,{g},"display","off");
    if numel(categories(g)) < 2
        warning('%s: need >=2 groups for Dunnett; skipping.',all_reg(i)); continue
    end
    [c,~]=multcompare(stats,"CriticalValueType","dunnett","ControlGroup",1);
    grpNames = string(stats.grpnames{1});
    for nc = 1:size(c,1)
        otherIdx = c(nc,2);
        if otherIdx == 1, otherIdx = c(nc,1); end
        idx = find(axonID(:)+" mu2" == grpNames(otherIdx), 1);
        if ~isempty(idx) && c(nc,6) > 0.05
            gofTable.boundaries{axonNum(idx)} = [];
        end
    end
end
%}
%% Per-axon bimodality summary
mu2Valid = gofTable.mu2Valid;

fprintf('\n--- mu2 validity breakdown (n = %d) ---\n',height(gofTable));
fprintf('  Fitted (K=2)         : %d\n',sum(gofTable.fitStatus=="fitted"));
fprintf('  mu2 > %.1f uV        : %d\n',P.minMu2_uV,sum(gofTable.thresholdPass));
fprintf('  Has boundary         : %d\n',sum(gofTable.hasBoundary));
fprintf('  dBIC > %g            : %d\n',P.minDeltaBIC,sum(gofTable.bicPass));
fprintf('  Weight/count pass    : %d\n',sum(gofTable.wPass));
fprintf('  Separation pass      : %d\n',sum(gofTable.sepPass));
fprintf('  ALL criteria         : %d\n',sum(mu2Valid));
for i = 1:length(all_reg)
    sel = gofTable.Subregion==all_reg(i) & hasLFP;
    fprintf('  %-10s %d/%d valid\n',all_reg(i),sum(sel & mu2Valid),sum(sel));
end

% Which criterion combinations are causing failures?
failSel = hasLFP & ~mu2Valid;
[reason,~,k] = unique(gofTable.failedCriteria(failSel));
disp(table(reason,accumarray(k,1),'VariableNames',{'FailedCriteria','Count'}))
%% More Statistical tests and figure generation

if ~exist("gofTable","var")
    gofTable = readtable( ...
        fullfile(figOutDir, ...
        "byAxonGMMFits_muSort_final.xlsx"));
end

%ampRange=[min(cell2mat(LFPTable.LFP),[],"all"),max(cell2mat(LFPTable.LFP),[],"all")];
ampRange = globalAmpRange;
binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
binCenters=convert_edges_2_centers(log10(binEdges))';
x_grid = linspace(log10(ampRange(1)),log10(ampRange(2)),100)';
gauss = @(x, m, s) (1 / sqrt(2*pi*s)) .* exp(-((x - m).^2) / (2*s));

%{
% ===== Rebuild per-axon validity from the Dunnett-pruned boundaries ======
% The Dunnett cell sets gofTable.boundaries{r} = [] for axons it rejects.
% An axon is therefore "valid" if it still has a boundary AND its mu2 sits
% above that boundary -- identical to the original iscell/mu2>boundaries
% logic, just collected once instead of branching inside every loop.
% Handles both the cell column (live gofTable) and the numeric column
% (gofTable restored via readtable).
mu2Valid = false(height(gofTable),1);
for r = 1:height(gofTable)
    if iscell(gofTable.boundaries)
        b = gofTable.boundaries{r};
    else
        b = gofTable.boundaries(r);
    end
    if ~isempty(b) && all(isfinite(b)) && gofTable.mu2(r) > b(1)
        mu2Valid(r) = true;
    end
end
fprintf('\nPost-Dunnett valid mu2: %d/%d axons\n', sum(mu2Valid), height(gofTable));
for i = 1:length(all_reg)
    sel = gofTable.Subregion == all_reg(i);
    fprintf('  %-10s %d/%d\n', all_reg(i), sum(sel & mu2Valid), sum(sel));
end
% ========================================================================
%}
% NEW: use the non-destructive per-electrode validity result.
mu2Valid = logical(gofTable.mu2Valid);

fprintf('\nValid mu2: %d/%d electrode pairs\n', ...
    sum(mu2Valid),height(gofTable));

for i = 1:length(all_reg)
    sel = gofTable.Subregion == all_reg(i) & hasLFP;

    fprintf('  %-10s %d/%d\n', ...
        all_reg(i), ...
        sum(sel & mu2Valid), ...
        sum(sel));
end
% Anova variables
x=[];
g=[];
xmu2=[];
gmu2=[];

for i=1:length(all_reg)
    figure
    hold on
    regRows = find(LFPTable.Subregion == all_reg(i) & hasLFP);
    if isempty(regRows)
        warning('%s has no usable LFP data; skipping.',all_reg(i));
        continue
    end
    axonNames=string();
    noiseFitMat=[];
    goodMu2=[];
    weight2=[];
    % [maxAmp,maxAmpIdx]=cellfun(@max,LFPTable.LFP(regRows));
    % [~,maxAmpAxonIdx]=max(maxAmp);
    % Obtain one scalar maximum per electrode.
    maxAmp = nan(numel(regRows),1);

    for k = 1:numel(regRows)
        thisLFP = LFPTable.LFP{regRows(k)};
        thisLFP = thisLFP(isfinite(thisLFP));

        if ~isempty(thisLFP)
            maxAmp(k) = max(thisLFP(:));
        end
    end

    if any(isfinite(maxAmp))
        [~,maxAmpAxonIdx] = max(maxAmp,[],'omitnan');
    else
        maxAmpAxonIdx = [];
    end

    % single loop -- each axon counted exactly once
    for naxons=1:length(regRows)
        r = regRows(naxons);

        if mu2Valid(r)
            goodMu2=[goodMu2,gofTable.mu2(r)];
            weight2=[weight2,gofTable.weights2(r)];
        end

        % Plot max amplitude axon
        %if naxons==maxAmpAxonIdx
        if ~isempty(maxAmpAxonIdx) && naxons == maxAmpAxonIdx
            histogram(log10(LFPTable.LFP{r}), log10(binEdges), 'Normalization','pdf', ...
                      'FaceColor',[0.2 0.6 1.0], 'FaceAlpha',0.5)
        end
    end

    % myMu1=gofTable.mu1(regRows);
    % mu1Avg = mean(myMu1);
    % newSD = std(gofTable.mu1(regRows));
    % meanPrior1=max(gofTable.weights1(regRows));
    validMu1Rows = regRows(isfinite(gofTable.muNoise(regRows)));
    myMu1 = gofTable.muNoise(validMu1Rows);
    if isempty(myMu1)
        warning('%s has no noise estimates.',all_reg(i)); close(gcf); continue
    end
    mu1Avg      = mean(myMu1,'omitnan');
    mu1Variance = var(myMu1,'omitnan');
    meanPrior1  = max(gofTable.wNoise(validMu1Rows),[],'omitnan');

    if isempty(myMu1)
        warning('%s has no successfully fitted mu1 values.',all_reg(i));
        close(gcf)
        continue
    end

    mu1Avg = mean(myMu1,'omitnan');
    mu1Variance = var(myMu1,'omitnan');
    meanPrior1 = max(gofTable.weights1(validMu1Rows),[],'omitnan');

    if isfinite(mu1Variance) && mu1Variance > 0 && isfinite(meanPrior1)
        plot(x_grid, ...
            gauss(x_grid,mu1Avg,mu1Variance).*meanPrior1, ...
            "LineWidth",4,"Color","b")
    end

    if ~isempty(goodMu2) && ~isempty(weight2)
        % mu2Avg = mean(goodMu2);
        % newSD = std(goodMu2);
        % if ~isfinite(newSD) || newSD == 0, newSD = eps; end   % single axon -> std = 0
        % meanPrior2 = max(weight2);
        % plot(x_grid,gauss(x_grid,mu2Avg,newSD).*meanPrior2,"LineWidth",4,"Color","r")
        % area(x_grid,gauss(x_grid,mu2Avg,newSD).*meanPrior2, ...
        %      'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none')
        mu2Avg = mean(goodMu2,'omitnan');
        newVariance = var(goodMu2,'omitnan');
        meanPrior2 = max(weight2,[],'omitnan');

        if ~isfinite(newVariance) || newVariance <= 0
            newVariance = eps;
        end

        mu2Curve = ...
            gauss(x_grid,mu2Avg,newVariance).*meanPrior2;

        plot(x_grid,mu2Curve, ...
            "LineWidth",4,"Color","r")

        area(x_grid,mu2Curve, ...
            'FaceColor','r', ...
            'FaceAlpha',0.25, ...
            'EdgeColor','none')
    else
        warning("%s: no valid mu2 to plot.", all_reg(i))
    end
    hold off

    ax=gca;
    xticks([-2:4])
    title(all_reg(i)+"  ("+sum(mu2Valid(regRows))+"/"+numel(regRows)+" valid mu2)")
    xlabel('Log_{10} \muV');
    ylabel('Probability Density');
    axis square
    box on
    set(ax,"FontSize",16)
    ax.LineWidth=4;
    ax.TickLength=[0.05,0.025];
    xlim(log10(ampRange))

    saveas(gcf, fullfile(figOutDir, "RegionPDF " + all_reg(i) + ".png"))

    x=[x;myMu1];
    g=[g;repmat(all_reg(i)+" mu1",[length(myMu1),1])];
    if ~isempty(goodMu2)
        x=[x;goodMu2'];
        g=[g;repmat(all_reg(i)+" mu2",[length(goodMu2),1])];
        xmu2=[xmu2;goodMu2'];
        gmu2=[gmu2;repmat(all_reg(i)+" mu2",[length(goodMu2),1])];
    end
end

% guarded ANOVAs -- a region with no valid mu2 contributes no mu2 group
if numel(unique(g)) >= 2
    figure
    [p,tbl,stats]=anovan(x,{g});
    [c,m]=multcompare(stats);
else
    warning('Fewer than 2 groups in g; skipping ANOVA.');
    p=[]; tbl=[]; stats=[]; c=[]; m=[];
end

if numel(unique(gmu2)) >= 2
    figure
    [pmu2,tblmu2,statsmu2]=anovan(xmu2,{gmu2});
    [cmu2,mmu2]=multcompare(statsmu2);
else
    warning('Only %d region(s) with valid mu2; skipping mu2 ANOVA.',numel(unique(gmu2)));
    pmu2=[]; tblmu2=[]; statsmu2=[]; cmu2=[]; mmu2=[];
end

fprintf('\n--- groups entering the ANOVA ---\n');
disp(unique(g))
fprintf('regions with a mu2 group: %d of %d\n', numel(unique(gmu2)), numel(all_reg));
%% Make bar chart

if isempty(stats)
    error('stats is empty -- the main ANOVA was skipped. Nothing to plot.');
end

grpNames = string(stats.grpnames{1});
means = nan(length(all_reg),2);      % rows = regions, col 1 = mu1, col 2 = mu2
SEs   = nan(length(all_reg),2);
nPerGroup = zeros(length(all_reg),2);

% OLD: do not reconstruct validity from boundaries alone.
%{
if ~exist('mu2Valid','var')
    mu2Valid = false(height(gofTable),1);
    for r = 1:height(gofTable)
        if iscell(gofTable.boundaries)
            b = gofTable.boundaries{r};
        else
            b = gofTable.boundaries(r);
        end
        if ~isempty(b) && all(isfinite(b)) && gofTable.mu2(r) > b(1)
            mu2Valid(r) = true;
        end
    end
end
%}

% NEW: use the same validity definition as the preceding analysis.
if ~exist('mu2Valid','var')
    if ismember("mu2Valid", ...
            string(gofTable.Properties.VariableNames))
        mu2Valid = logical(gofTable.mu2Valid);
    else
        error(['gofTable.mu2Valid is missing. Run the per-axon ' ...
            'bimodality-criteria section before this section.']);
    end
end

for i = 1:length(all_reg)
    i1 = find(grpNames == all_reg(i)+" mu1",1);
    i2 = find(grpNames == all_reg(i)+" mu2",1);
    if ~isempty(i1), means(i,1)=m(i1,1); SEs(i,1)=m(i1,2); end
    if ~isempty(i2), means(i,2)=m(i2,1); SEs(i,2)=m(i2,2); end

    sel = gofTable.Subregion == all_reg(i);

    % OLD: this counted skipped and failed rows.
    % nPerGroup(i,1) = sum(sel);

    % NEW: count only the finite fitted mu1 values that entered the ANOVA.
    validMu1Selection = sel & isfinite(gofTable.muNoise);

    nPerGroup(i,1) = sum(validMu1Selection);
    nPerGroup(i,2) = sum(sel & mu2Valid);
end

myCats=categorical(all_reg);
myCats=reordercats(myCats,all_reg);

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
set(gcf, 'Toolbar', 'none', 'Menu', 'none')
hBar=bar(myCats, means);
hold on

nSeries = size(means,2);             % 2 (mu1, mu2)
for s = 1:nSeries
    xe = hBar(s).XEndPoints;
    for gi = 1:numel(xe)
        if isfinite(means(gi,s))
            errorbar(xe(gi), means(gi,s), SEs(gi,s), 'k', 'linestyle','none', 'CapSize',10);
        end
    end
end

% Add sample counts above the error bars.
yl_now = ylim;
labelOffset = 0.03*range(yl_now);

for s = 1:nSeries
    xe = hBar(s).XEndPoints;

    for gi = 1:numel(xe)
        if isfinite(means(gi,s))
            text(xe(gi), ...
                means(gi,s) + SEs(gi,s) + labelOffset, ...
                string(nPerGroup(gi,s)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',12);
        else
            text(xe(gi), ...
                yl_now(1) + labelOffset, ...
                "0", ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',12, ...
                'Color',[0.5 0.5 0.5]);
        end
    end
end

% Add vertical headroom so the counts are not clipped.
yl = ylim;
ylim([yl(1), yl(2) + 0.12*range(yl)])

hold off

title("Low Gamma Source Amplitudes Above Noise")
%title("High Gamma Source Amplitudes Above Noise")
subtitle("Numbers above bars indicate n (fitted electrode–axon pairs)")
xlabel("Subregion")
ylabel("Average log_{10} \muV")
legend({'\mu_1 (noise)','\mu_2 (signal)'},'Location','eastoutside')

ax=gca;
ax.FontSize=16;
axis square
box on
ax.LineWidth=4;
ax.TickLength=[0.05,0.025];

saveas(gcf, fullfile(figOutDir, "Mu1Mu2_BarChart.png"))

%% ===================== Local functions =====================
function [best,src] = fitBestTwoComponent(x,P,options,regValue)
best = []; src = "";
cands = {}; names = strings(0);
lowD  = x(x <= P.minMu2_log);
highD = x(x >  P.minMu2_log);

% Targeted start only when both sides of the threshold have enough data
if numel(highD) >= P.minCount2 && numel(lowD) >= P.minCount2
    s = struct();
    s.mu = [median(lowD); median(highD)];
    v1 = max(var(lowD,1),1e-4);
    v2 = max(var(highD,1),(1.5^2)*v1);
    s.Sigma = reshape([v1 v2],1,1,2);
    w2 = min(max(mean(x > P.minMu2_log),0.05),0.40);
    s.ComponentProportion = [1-w2, w2];
    try
        cands{end+1} = fitgmdist(x,2,'Start',s,'Options',options, ...
            'RegularizationValue',regValue);
        names(end+1) = "targeted";
    catch
    end
end

% Automatic k-means++ restarts (always)
try
    cands{end+1} = fitgmdist(x,2,'Replicates',P.nReplicates, ...
        'Options',options,'RegularizationValue',regValue);
    names(end+1) = "kmeans++";
catch
end

if isempty(cands), return; end
nll = cellfun(@(g) g.NegativeLogLikelihood,cands);
[~,k] = min(nll);
best = cands{k}; src = names(k);
end

function bnd = findBetweenMeanBoundary(mu,vr,w,xMin,xMax)
bnd = NaN;
if any(~isfinite([mu(:);vr(:);w(:)])), return; end
x = linspace(min(xMin,mu(1)),max(xMax,mu(2)),20000)';
d = w(1)*normpdf(x,mu(1),sqrt(vr(1))) - w(2)*normpdf(x,mu(2),sqrt(vr(2)));
k = find(d(1:end-1).*d(2:end) <= 0);
if isempty(k), return; end
xc = x(k) - d(k).*(x(k+1)-x(k))./(d(k+1)-d(k));
bad = ~isfinite(xc); xc(bad) = x(k(bad));
xc = xc(xc > mu(1) & xc < mu(2));
if isempty(xc), return; end
[~,j] = min(abs(xc-mean(mu)));
bnd = xc(j);
end

function C = evaluateMu2Criteria(mu,sd,w,bnd,bic1,bic2,N,P)
C.threshold  = isfinite(mu(2)) && mu(2) > P.minMu2_log;
C.boundary   = isfinite(bnd) && bnd > mu(1) && bnd < mu(2);
C.bic        = isfinite(bic1) && isfinite(bic2) && (bic1-bic2) > P.minDeltaBIC;
C.weight     = isfinite(w(2)) && w(2) >= P.minWeight2 && w(2)*N >= P.minCount2;
C.separation = isfinite(sd(1)) && sd(1) > 0 && mu(2) > mu(1) + P.minSeparation*sd(1);
boundaryOK   = C.boundary || ~P.requireBoundary;
C.valid = C.threshold && boundaryOK && C.bic && C.weight && C.separation;

f = strings(0,1);
if ~C.threshold,             f(end+1) = "mu2 below " + P.minMu2_uV + " uV"; end
if ~C.boundary && P.requireBoundary, f(end+1) = "no boundary"; end
if ~C.bic,                   f(end+1) = "BIC"; end
if ~C.weight,                f(end+1) = "weight/count"; end
if ~C.separation,            f(end+1) = "separation"; end
if isempty(f), C.failedText = "none"; else, C.failedText = strjoin(f,", "); end
end

function fig = plotGmmDiagnostic(x,edgesLog,name,mu,sd,w,bnd,C,muK1,sdK1,P)
fig = figure('Visible','off','Color','w');
hold on
histogram(x,edgesLog,'Normalization','pdf', ...
    'FaceColor',[0.55 0.78 1.00],'EdgeColor',[0.10 0.25 0.40], ...
    'FaceAlpha',0.65,'DisplayName','Observed data');
pad = max(0.1*range(x),0.1);
xg = linspace(min(x)-pad,max(x)+pad,2000)';

if C.valid
    j1 = w(1)*normpdf(xg,mu(1),sd(1));
    j2 = w(2)*normpdf(xg,mu(2),sd(2));
    plot(xg,j1,'b-','LineWidth',3,'DisplayName',sprintf('\\mu_1 (noise) = %.3f',mu(1)));
    plot(xg,j2,'r-','LineWidth',3,'DisplayName',sprintf('\\mu_2 (signal) = %.3f',mu(2)));
    plot(xg,j1+j2,'k--','LineWidth',2,'DisplayName','2-component mixture');
    xline(mu(1),'b:','LineWidth',2,'HandleVisibility','off');
    xline(mu(2),'r:','LineWidth',2,'HandleVisibility','off');
    if isfinite(bnd)
        xline(bnd,'m-','LineWidth',3,'DisplayName',sprintf('Boundary = %.3f',bnd));
    end
    ttl = sprintf('PASS | \\mu_1=%.3f, \\mu_2=%.3f, SD_1=%.3f, SD_2=%.3f, w=%.2f/%.2f', ...
        mu(1),mu(2),sd(1),sd(2),w(1),w(2));
else
    if P.showRejectedMixture && all(isfinite([mu;sd;w]))
        plot(xg,w(1)*normpdf(xg,mu(1),sd(1))+w(2)*normpdf(xg,mu(2),sd(2)), ...
            '-','Color',[0.5 0.5 0.5],'LineWidth',1,'DisplayName','Rejected 2-comp mixture');
    end
    plot(xg,normpdf(xg,muK1,sdK1),'b-','LineWidth',3, ...
        'DisplayName',sprintf('\\mu_1 (single Gaussian) = %.3f',muK1));
    xline(muK1,'b:','LineWidth',2,'HandleVisibility','off');
    ttl = sprintf('FAIL | \\mu_1=%.3f, SD_1=%.3f | failed: %s', ...
        muK1,sdK1,char(C.failedText));
end

title({char(name),ttl});
xlim(edgesLog([1 end]))
xlabel('Log_{10} \muV'); ylabel('Probability Density');
legend('show','Location','best'); box on; hold off
end