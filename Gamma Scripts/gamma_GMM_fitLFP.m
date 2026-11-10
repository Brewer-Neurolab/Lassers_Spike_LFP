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

    xline(ax1,m+2*v)
    legend(ax1,axonNames)
    xline(ax2,m+2*v)

    % ===== EDIT 4: clamp last_bin. If m+2*v lands past the final edge (common
    %               in high gamma, short tail) last_bin(1) errored on empty.
    % last_bin=find(binEdges>m+2*v);
    % last_bin=last_bin(1);
    last_bin=find(binEdges>m+2*v);
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
    xline(ax3,m+2*v,'LineWidth',4,'Color','k','LineStyle','--')
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
    xline(ax4,m+2*v,'LineWidth',4,'Color','k','LineStyle','--')
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

%% Gaussian Mixed Model
close all
rng('default')
gofTable=LFPTable(:,1:6);
% Explicitly initialize cell columns for every table row.
gofTable.gmm = cell(height(gofTable),1);
gofTable.myData = cell(height(gofTable),1);
gofTable.boundaries = cell(height(gofTable),1);
gofTable.originalBounds = cell(height(gofTable),1);
% Track processing status without changing the fit logic.
gofTable.fitStatus = repmat("not_processed",height(gofTable),1);
gofTable.skipReason = repmat("",height(gofTable),1);

% Initialize t-test outputs as missing.
% NaN means that the test was not performed.
gofTable.ttest_h = nan(height(gofTable),1);
gofTable.ttest_p = nan(height(gofTable),1);
gofTable.ttestPerformed = false(height(gofTable),1);

% Optional diagnostic for extremely small p-values.
gofTable.ttest_negLog10P = nan(height(gofTable),1);

gofRow=1;

nClasses=2;

% ampRange=[min(cell2mat(LFPTable.LFP(gofRow)),[],"all"),max(cell2mat(LFPTable.LFP(gofRow)),[],"all")];
% binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
% binCenters=convert_edges_2_centers(log10(binEdges))';
% ampRange=[min(cell2mat(LFPTable.LFP),[],"all"),max(cell2mat(LFPTable.LFP),[],"all")];
% binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
% binCenters=convert_edges_2_centers(log10(binEdges))';
% x_grid = linspace(log10(ampRange(1)),log10(ampRange(2)),10000)';
% Use the finite positive amplitudes collected from usable rows.
ampRange = globalAmpRange;
binEdges = logspace(log10(ampRange(1)),log10(ampRange(2)),101);
binCenters = convert_edges_2_centers(log10(binEdges))';
x_grid = linspace(log10(ampRange(1)),log10(ampRange(2)),10000)';

% Fit every row that contains usable LFP amplitudes.
for i = 1:numel(processedRows)

    % Store the result in its original LFPTable/gofTable row.
    gofRow = processedRows(i);

    % regRows=find(gofTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    
    ampPts=[];
    axonNames=string();
    % for naxons=1:length(regRows)
    figure
    hold on
    myData = log10(LFPTable.LFP{gofRow});
    myData = myData(:);   % force column vector regardless of input orientation
    % ===== EDIT 6: drop any residual -Inf/NaN before fitting.
    myData = myData(isfinite(myData));
    % ======================================================================
    pts=histogram(myData,log10(binEdges));
    ampPts=pts.BinCounts;
    axonNames(i)=gofTable.Subregion(gofRow)+" FID "+gofTable.fi(gofRow)+" "+gofTable.Electrode(gofRow);
    hold off

    % % Fit a 2-Component GMM Using EM
    % % mu = [1;100];
    % % Sigma(:,:,1) = 1;
    % % Sigma(:,:,2) = 2;
    % % % start = struct('mu', mu, 'Sigma', Sigma);
    % options = statset('MaxIter', 1000, 'TolFun', 1e-8);
    % start.mu = [quantile(myData, 0.25); quantile(myData, 0.75)];
    % start.Sigma = cat(3, var(myData)*0.5, var(myData)*0.5);
    % gmm = fitgmdist(myData, nClasses, 'Options', options, ...
    %             'RegularizationValue', 1e-4, ...
    %             'Replicates', 1, ...
    %             'Start', start);
    % % % --- NEW: also fit K=1 so we can ask whether 2 components are justified ---
    % gmm1 = fitgmdist(myData, 1, 'Options', options, 'Replicates', 1, ...
    %                  'RegularizationValue', 1e-6);
    % % ------------------------------------------------------------------------
    % % 'Replicates' runs EM from 5 different starting points
    % % to reduce the chance of landing on a local minimum
    % 
    % % Extract GMM Parameters
    % mu      = gmm.mu;                        % 2x1 means
    % sigma   = squeeze(gmm.Sigma);            % 2x1 variances (squeezed from 1x1x2)
    % weights = gmm.ComponentProportion;       % 1x2 mixture weights
    % 
    % % Calculate amplitudes
    % amp_unweighted=[];
    % amp_weighted=[];
    % for nAmps = 1:nClasses
    %     amp_unweighted(nAmps) = 1 / (sqrt(2 * pi) * sqrt(sigma(nAmps)));
    %     amp_weighted(nAmps)   = weights(nAmps) * amp_unweighted(nAmps);
    %     fprintf('Component %d: Amplitude = %.4f (unweighted), %.4f (weighted)\n', ...
    %         nAmps, amp_unweighted(nAmps), amp_weighted(nAmps));
    % end
    % 
    % % ===== OLD (comment out this whole block) =============================
    % % % Sort amplitudes in descending order
    % % [weightedAmp_sorted, idx] = sort(amp_weighted,"descend");
    % % unweightedAmp_sorted = amp_unweighted(idx);
    % % mu_sorted = mu(idx);
    % % sigma_sorted  = sigma(idx);
    % % weights_sorted = weights(idx);
    % %
    % % % Assign sorted parameters
    % % amp_weighted = weightedAmp_sorted;
    % % amp_unweighted = unweightedAmp_sorted;
    % % mu      = mu_sorted;
    % % sigma   = sigma_sorted;
    % % weights = weights_sorted;
    % % ======================================================================
    % 
    % % ===== NEW: sort by mu ASCENDING ======================================
    % % amp_weighted = w/(sqrt(2*pi)*sigma) is a HEIGHT, not a POSITION, so a
    % % narrow low-weight component can outrank a broad high-weight one and the
    % % labels invert. In byAxonGMMFits.xlsx roughly half the rows are inverted,
    % % e.g. EC-DG FID 4 F2: mu1 = 0.19610, mu2 = -0.09938 (comp 2 BELOW comp 1).
    % % Sorting by mu makes mu1 the low/noise mode and mu2 the high/signal mode
    % % on every axon, by definition.
    % [mu, idx]      = sort(mu, 'ascend');
    % sigma          = sigma(idx);
    % weights        = weights(idx);
    % amp_unweighted = amp_unweighted(idx);
    % amp_weighted   = amp_weighted(idx);
    % % ======================================================================
    % 
    % fprintf('Component 1: mu = %.4f, sigma^2 = %.4f, weight = %.4f\n', ...
    %     mu(1), sigma(1), weights(1));
    % fprintf('Component 2: mu = %.4f, sigma^2 = %.4f, weight = %.4f\n', ...
    %     mu(2), sigma(2), weights(2));
    % 
    % % Class Priors
    % prior1 = weights(1);
    % prior2 = weights(2);
    % 
    % % Define Likelihoods and Posteriors
    % % Gaussian PDF for each component
    % gauss = @(x, m, s) (1 / sqrt(2*pi*s)) .* exp(-((x - m).^2) / (2*s));
    % 
    % % Evaluate over a fine grid for plotting and boundary finding
    % % x_grid = linspace(min(myData)-0.5, max(myData)+0.5, 10000)';
    % 
    % % Likelihoods: P(x | class)
    % likelihood1 = gauss(x_grid, mu(1), sigma(1));
    % likelihood2 = gauss(x_grid, mu(2), sigma(2));
    % 
    % % Joint probabilities: P(x | class) * P(class)
    % joint1 = prior1 .* likelihood1;
    % joint2 = prior2 .* likelihood2;
    % 
    % % Evidence: P(x) — used to normalize posteriors
    % evidence = joint1 + joint2;
    % 
    % % Posteriors: P(class | x) via Bayes' theorem
    % posterior1 = joint1 ./ evidence;
    % posterior2 = joint2 ./ evidence;
    % 
    % % Find the Bayesian Decision Boundary
    % % The boundary is where posterior1 == posterior2, i.e., difference crosses 0
    % diff_posteriors = posterior1 - posterior2;
    % sign_changes = find(diff(sign(diff_posteriors)) ~= 0);
    % 
    % % There may be one or two crossing points depending on overlap
    % % CONSIDER ONLY POSITIVE BOUNDARIES AND BOUNDARIES BETWEEN BOTH MU
    % boundaries = x_grid(sign_changes);
    % originalBounds=boundaries;
    % % boundaries = boundaries(boundaries>0);
    % maxAmpIdx=find(amp_weighted==max(amp_weighted));
    % 
    % % boundaryThresh=mu(maxAmpIdx(1))+sqrt(sigma(maxAmpIdx(1)));
    % % boundaries = boundaries(boundaries>boundaryThresh);
    % % if length(boundaries)>1
    % %     dBound1=abs(sum(boundaries(1)-mu));
    % %     dBound2=abs(sum(boundaries(2)-mu));
    % %     [~,dBoundMinIdx]=min([dBound1,dBound2]);
    % %     boundaries=boundaries(dBoundMinIdx);
    % % end
    % % boundaries = boundaries(boundaries<=max(mu) & boundaries>=min(mu));
    % 
    % % ===== OLD (comment out) ==============================================
    % % % mu is already sorted descending by weighted amplitude, so index 1 = dominant/noise
    % % boundaryThresh = mu(1);                                  % was mu(1)+1*SD -- too strict for gamma
    % % boundaries = boundaries(boundaries>boundaryThresh & boundaries<max(myData));
    % % if isempty(boundaries)                                   % NEW: single-Gaussian-like fallback
    % %     boundaries = mu(1) + 3*sqrt(sigma(1));               % 3 SD above the noise mode
    % % end
    % % boundaries = min(boundaries);                            % force scalar
    % % ======================================================================
    % 
    % % ===== NEW: boundary strictly between the two means, never fabricated ==
    % % mu is now sorted ASCENDING, so mu(1) = low mode, mu(2) = high mode.
    % % Two components always cross twice; only the crossing BETWEEN the means
    % % is a meaningful low-vs-high threshold. The other crossing is just the
    % % broader component's tail overtaking the narrower one, which is what
    % % produced x = 0.695 on EC-DG FID 4 F2 -- a line sitting in empty space
    % % past the end of that histogram.
    % % ===== [mu sort block ends here] =====
    % 
    % % Fix 2: minimum separation check
    % muSeparation = mu(2) - mu(1);
    % minSeparation = 0.5;   % log10 units
    % 
    % if muSeparation < minSeparation
    %     boundaries = [];
    %     disp("WARNING: mu separation too small (" + muSeparation + ") for " + axonNames(i))
    % end
    % 
    % boundaryThresh = mu(1);
    % 
    % % Stage 1: ideal case — boundary strictly between the two means
    % boundaries = originalBounds(originalBounds > mu(1) & originalBounds < mu(2));
    % 
    % if ~isempty(boundaries)
    %     boundaries = boundaries(1);   % leftmost crossing between modes
    % else
    %     % Stage 2: fallback for overlapping/skewed distributions —
    %     % take the leftmost crossing to the right of mu(1) up to max(myData).
    %     % This recovers the tail-crossing case (e.g. CA1-EC FID 1 A7).
    %     fallback = originalBounds(originalBounds > mu(1) & originalBounds < max(myData));
    %     if ~isempty(fallback)
    %         boundaries = min(fallback);
    %     else
    %         boundaries = [];   % genuinely no separating boundary
    %     end
    % end
    % % ======================================================================
    % 
    % disp("Min Boundary X Location: "+boundaryThresh)
    % 
    % fprintf('\nBayesian Decision Boundary/Boundaries:\n');
    % for nBoundaries = 1:length(boundaries)
    %     fprintf('  Boundary %d: x = %.4f\n', nBoundaries, boundaries(nBoundaries));
    % end
    % 
    % % Classify Data Using the Boundary
    % % For a single boundary, classify based on side
    % if length(boundaries) == 1
    %     predicted_labels = double(myData > boundaries(1)) + 1;
    %     % Class 1 = left of boundary, Class 2 = right of boundary
    % elseif ~isempty(boundaries)
    %     % For two boundaries, assign based on nearest component mean
    %     predicted_labels = ones(size(myData));
    %     for nData = 1:length(myData)
    %         p1 = prior1 * gauss(myData(nData), mu(1), sigma(1));
    %         p2 = prior2 * gauss(myData(nData), mu(2), sigma(2));
    %         if p2 > p1
    %             predicted_labels(nData) = 2;
    %         end
    %     end
    % end
    
    % -------------------------------------------------------------------------
    % Fit a two-component GMM using multiple initializations
    % -------------------------------------------------------------------------
    options = statset('MaxIter',1000,'TolFun',1e-8);

    % OLD: fixed start with only one replicate
    % start.mu = [quantile(myData,0.25); quantile(myData,0.75)];
    % start.Sigma = cat(3,var(myData)*0.5,var(myData)*0.5);
    % gmm = fitgmdist(myData,nClasses,'Options',options, ...
    %     'RegularizationValue',1e-4, ...
    %     'Replicates',1, ...
    %     'Start',start);

    % Small scale-dependent regularization prevents collapsed variances.
    dataVariance = var(myData,1);
    regValue = max(1e-10,1e-6*dataVariance);

    % % Fifty automatic initializations; no custom starting structure.
    % gmm = fitgmdist(myData,2, ...
    %     'Options',options, ...
    %     'Replicates',50, ...
    %     'RegularizationValue',regValue);

    % -------------------------------------------------------------------------
    % Targeted initialization:
    % Start component 2 above 1 on the log10 plot (>10 uV), and initialize
    % component 2 with a larger standard deviation than component 1.
    % -------------------------------------------------------------------------
    minMu2_log = log10(30);   % log10 amplitude of 1 = 10 uV

    lowStartData  = myData(myData <= minMu2_log);
    highStartData = myData(myData >  minMu2_log);

    % -------------------------------------------------------------------------
    % Skip electrodes that cannot support the requested initialization.
    % -------------------------------------------------------------------------
    if isempty(highStartData) || isempty(lowStartData)

        maxAmplitude_uV = 10^max(myData);

        if isempty(highStartData)
            skipReason = "no_data_above_10uV";

            warning(['Skipping two-component GMM for %s: no data exceed ' ...
                '10 uV. Maximum amplitude = %.3f uV.'], ...
                axonNames(i),maxAmplitude_uV);
        else
            skipReason = "no_data_at_or_below_10uV";

            warning(['Skipping two-component GMM for %s: no data are at or ' ...
                'below 10 uV.'],axonNames(i));
        end

        % Record why this row was not fitted.
        gofTable.fitStatus(gofRow) = "skipped";
        gofTable.skipReason(gofRow) = skipReason;
        gofTable.nSamples(gofRow) = numel(myData);
        gofTable.maxAmplitude_uV(gofRow) = maxAmplitude_uV;

        % Explicit missing values prevent skipped rows from appearing to have
        % valid zero-valued fits.
        gofTable.NegativeLogLikelihood(gofRow) = NaN;
        gofTable.BIC1(gofRow) = NaN;
        gofTable.BIC2(gofRow) = NaN;

        gofTable.mu1(gofRow) = NaN;
        gofTable.mu2(gofRow) = NaN;
        gofTable.sigma1(gofRow) = NaN;
        gofTable.sigma2(gofRow) = NaN;
        gofTable.std1(gofRow) = NaN;
        gofTable.std2(gofRow) = NaN;
        gofTable.weights1(gofRow) = NaN;
        gofTable.weights2(gofRow) = NaN;

        gofTable.integration_log10(gofRow) = NaN;
        gofTable.integration(gofRow) = NaN;
        gofTable.percentIntegrated(gofRow) = NaN;

        gofTable.boundaries{gofRow} = [];
        gofTable.originalBounds{gofRow} = [];
        gofTable.gmm{gofRow} = [];
        gofTable.myData{gofRow} = myData;

        % Close the temporary histogram figure created earlier in this iteration.
        if ~isempty(get(groot,'CurrentFigure'))
            close(gcf);
        end

        % Advance to the next table row and skip the remaining GMM code.
        %gofRow = gofRow+1;
        continue
    end

    % Clear the initialization from the previous electrode.
    clear start

    % Start mu1 in the main lower-amplitude distribution.
    start.mu(1,1) = median(lowStartData);

    % Start mu2 in the high-amplitude data, strictly above 1 log10 unit.
    start.mu(2,1) = median(highStartData);
    start.mu(2,1) = max(start.mu(2,1),minMu2_log + 0.01);

    % Estimate the lower-component starting variance.
    var1Start = var(lowStartData,1);

    if ~isfinite(var1Start) || var1Start <= 0
        var1Start = max(0.25*var(myData,1),1e-4);
    end

    % Estimate the high-component variance, but force its starting SD to be
    % larger than component 1. A 1.5-times larger SD means a 2.25-times
    % larger variance.
    var2Observed = var(highStartData,1);

    if ~isfinite(var2Observed) || var2Observed <= 0
        var2Observed = var1Start;
    end

    var2Start = max(var2Observed,(1.5^2)*var1Start);

    start.Sigma = zeros(1,1,2);
    start.Sigma(:,:,1) = var1Start;
    start.Sigma(:,:,2) = var2Start;

    % Initialize the high-amplitude component using the observed proportion
    % above 10 uV, while preventing an extremely small initial component.
    highStartWeight = mean(myData > minMu2_log);
    highStartWeight = min(max(highStartWeight,0.05),0.40);

    start.ComponentProportion = ...
        [1-highStartWeight, highStartWeight];

    % A custom Start structure specifies one targeted initialization.
    gmm = fitgmdist(myData,2, ...
        'Options',options, ...
        'Start',start, ...
        'Replicates',1, ...
        'RegularizationValue',regValue);

    % OLD: do not mark the row fitted until both models finish.
    % gofTable.fitStatus(gofRow) = "fitted";
    % gofTable.skipReason(gofRow) = "";

    % Keep this only if you want BIC comparison with a one-component model.
    gmm1 = fitgmdist(myData,1, ...
        'Options',options, ...
        'RegularizationValue',regValue);

    % Both requested fits completed successfully.
    gofTable.fitStatus(gofRow) = "fitted";
    gofTable.skipReason(gofRow) = "";

    % -------------------------------------------------------------------------
    % Sort components by mean
    % Component 1 = lower-amplitude mode
    % Component 2 = higher-amplitude mode
    % -------------------------------------------------------------------------
    [mu,idx] = sort(gmm.mu(:),'ascend');

    sigmaAll  = squeeze(gmm.Sigma);
    sigmaAll  = sigmaAll(:);   % These values are variances
    weightsAll = gmm.ComponentProportion(:);

    sigma  = sigmaAll(idx);
    weights = weightsAll(idx);
    
    % -------------------------------------------------------------------------
    % Check whether the final solution retained the targeted properties.
    % sigma contains variance, so take the square root to obtain SD.
    % -------------------------------------------------------------------------
    stdComponents = sqrt(sigma);

    finalMu2AboveOne = mu(2) > minMu2_log;
    finalStd2Larger  = stdComponents(2) > stdComponents(1);

    fprintf(['Start: mu1 = %.4f, mu2 = %.4f, ' ...
        'SD1 = %.4f, SD2 = %.4f\n'], ...
        start.mu(1),start.mu(2), ...
        sqrt(start.Sigma(:,:,1)),sqrt(start.Sigma(:,:,2)));

    fprintf(['Final: mu1 = %.4f, mu2 = %.4f, ' ...
        'SD1 = %.4f, SD2 = %.4f\n'], ...
        mu(1),mu(2),stdComponents(1),stdComponents(2));

    if ~finalMu2AboveOne
        warning('Final mu2 moved below 1 log10 unit for %s.',axonNames(i));
    end

    if ~finalStd2Larger
        warning('Final SD2 is not larger than SD1 for %s.',axonNames(i));
    end

    % Peak heights after sorting
    amp_unweighted = 1./sqrt(2*pi*sigma);
    amp_weighted   = weights.*amp_unweighted;

    fprintf('Component 1: mu = %.4f, variance = %.4f, weight = %.4f\n', ...
        mu(1),sigma(1),weights(1));

    fprintf('Component 2: mu = %.4f, variance = %.4f, weight = %.4f\n', ...
        mu(2),sigma(2),weights(2));

    % -------------------------------------------------------------------------
    % Evaluate the two component curves
    % -------------------------------------------------------------------------
    xPad = 0.10*range(myData);

    if xPad <= 0 || ~isfinite(xPad)
        xPad = 0.1;
    end

    x_grid = linspace(min(myData)-xPad,max(myData)+xPad,10000)';

    gauss = @(x,m,v) exp(-0.5*((x-m).^2)./v)./sqrt(2*pi*v);

    % Weighted component densities
    joint1 = weights(1).*gauss(x_grid,mu(1),sigma(1));
    joint2 = weights(2).*gauss(x_grid,mu(2),sigma(2));

    % Complete mixture density
    evidence = joint1+joint2;

    % Retain these variables for compatibility with your later code
    likelihood1 = gauss(x_grid,mu(1),sigma(1));
    likelihood2 = gauss(x_grid,mu(2),sigma(2));

    posterior1 = joint1./max(evidence,realmin);
    posterior2 = joint2./max(evidence,realmin);

    % -------------------------------------------------------------------------
    % Find only the decision boundary between the two means
    % -------------------------------------------------------------------------
    densityDifference = joint1-joint2;

    crossingIdx = find( ...
        densityDifference(1:end-1).*densityDifference(2:end) <= 0);

    originalBounds = [];

    % Interpolate each crossing for greater precision
    for kCross = 1:numel(crossingIdx)
        k = crossingIdx(kCross);

        x1 = x_grid(k);
        x2 = x_grid(k+1);
        y1 = densityDifference(k);
        y2 = densityDifference(k+1);

        if y2 ~= y1
            crossing = x1-y1*(x2-x1)/(y2-y1);
        else
            crossing = mean([x1,x2]);
        end

        originalBounds(end+1,1) = crossing;
    end

    % A meaningful low-vs-high boundary must fall between the two means.
    boundaries = originalBounds( ...
        originalBounds > mu(1) & originalBounds < mu(2));

    if ~isempty(boundaries)
        % Select the crossing nearest the midpoint of the two means.
        [~,nearestIdx] = min(abs(boundaries-mean(mu)));
        boundaries = boundaries(nearestIdx);
    else
        boundaries = [];
        warning('No decision boundary exists between mu1 and mu2.');
    end

    % Retained for compatibility with your table code
    boundaryThresh = mu(1);

    % Classify observations only if a between-mean boundary was found
    if ~isempty(boundaries)
        predicted_labels = double(myData > boundaries(1))+1;
    else
        predicted_labels = nan(size(myData));
    end

    % -------------------------------------------------------------------------
    % Perform the t-test only when a boundary produced two usable groups.
    %
    % IMPORTANT: this test is descriptive only. The same fitted boundary is
    % used to create the groups, so the test is not independent evidence of
    % bimodality.
    % -------------------------------------------------------------------------
    h = NaN;
    p = NaN;
    ttestPerformed = false;

    if ~isempty(boundaries)

        data1 = myData(predicted_labels == 1);
        data2 = myData(predicted_labels == 2);

        if numel(data1) > 1 && numel(data2) > 1
            [h,p] = ttest2(data1,data2, ...
                'Vartype','unequal', ...
                'Alpha',0.05);

            ttestPerformed = true;
        else
            warning(['The decision boundary did not produce two groups with ' ...
                'at least two observations for %s.'],axonNames(i));
        end
    end

    % Save the result for every row, including NaN when no test was performed.
    gofTable.ttest_h(gofRow) = h;
    gofTable.ttest_p(gofRow) = p;
    gofTable.ttestPerformed(gofRow) = ttestPerformed;

    % Store -log10(p), which remains interpretable for very small p-values.
    % If ttest2 returns exactly zero because of numerical underflow, this is Inf.
    if ttestPerformed
        if p > 0
            gofTable.ttest_negLog10P(gofRow) = -log10(p);
        elseif p == 0
            gofTable.ttest_negLog10P(gofRow) = Inf;
        end
    end

    % Figure rule:
    % Show mu2 only when the test was performed and rejected H0.
    % With no boundary, showMu2 is false because no test was performed.
    showMu2 = ttestPerformed && (h == 1);
    figure;
    hold on;

    % % Histogram of raw data
    % histogram(myData, log10(binEdges), 'Normalization', 'pdf', ...
    %     'FaceColor', [0.2 0.6 1.0], 'FaceAlpha', 0.5, 'DisplayName', 'Class 1 Data');
    % 
    % % GMM component PDFs (scaled by prior)
    % plot(x_grid, joint1, 'b-', 'LineWidth', 2, 'DisplayName', 'GMM Component 1');
    % plot(x_grid, joint2, 'r-', 'LineWidth', 2, 'DisplayName', 'GMM Component 2');
    % 
    % % Overall GMM mixture PDF
    % plot(x_grid, evidence, 'k--', 'LineWidth', 1.5, 'DisplayName', 'GMM Mixture');
    % 
    % % Decision boundary line(s)
    % for nBounds = 1:length(boundaries)
    %     xline(boundaries(nBounds), 'm-', 'LineWidth', 2.5, ...
    %         'DisplayName', sprintf('Decision Boundary: x=%.3f', boundaries(nBounds)));
    % end
    % 
    % xlabel('Log_{10} \muV');
    % ylabel('Probability Density');
    % title(axonNames(i));
    % legend('show', 'Location', 'best');
    % hold off;

    % -------------------------------------------------------------------------
    % Plot histogram, sorted components, mixture, means, and boundary
    % -------------------------------------------------------------------------
    figure;
    hold on;

    histogram(myData,log10(binEdges), ...
        'Normalization','pdf', ...
        'FaceColor',[0.55 0.78 1.00], ...
        'EdgeColor',[0.10 0.25 0.40], ...
        'FaceAlpha',0.65, ...
        'DisplayName','Observed Data');

    % Lower-mean component is always blue
    plot(x_grid,joint1,'b-', ...
        'LineWidth',3, ...
        'DisplayName',sprintf('Low-\\mu Component: \\mu=%.3f',mu(1)));

    % Higher-mean component is always red
    plot(x_grid,joint2,'r-', ...
        'LineWidth',3, ...
        'DisplayName',sprintf('High-\\mu Component: \\mu=%.3f',mu(2)));

    % Complete fitted distribution
    plot(x_grid,evidence,'k--', ...
        'LineWidth',2.5, ...
        'DisplayName','GMM Mixture');

    % Show the location of each component mean
    xline(mu(1),'b:', ...
        'LineWidth',2, ...
        'HandleVisibility','off');

    xline(mu(2),'r:', ...
        'LineWidth',2, ...
        'HandleVisibility','off');

    % Plot a boundary only when it occurs between the two means
    if ~isempty(boundaries)
        xline(boundaries(1),'m-', ...
            'LineWidth',3, ...
            'DisplayName',sprintf('Decision Boundary: x=%.3f',boundaries(1)));
    end

    xlabel('Log_{10} \muV');
    ylabel('Probability Density');

    title({axonNames(i), ...
        sprintf(['\\mu_1=%.3f, \\mu_2=%.3f, ' ...
        'SD_1=%.3f, SD_2=%.3f, weights=%.2f/%.2f'], ...
        mu(1),mu(2),stdComponents(1),stdComponents(2), ...
        weights(1),weights(2))});

    legend('show','Location','best');
    box on;
    hold off;

    gofTable.gmm{gofRow}=gmm;
    gofTable.NegativeLogLikelihood(gofRow)=gmm.NegativeLogLikelihood;
    % --- NEW: model-comparison bookkeeping ----------------------------------
    gofTable.BIC1(gofRow)     = gmm1.BIC;        % 1-component BIC
    gofTable.BIC2(gofRow)     = gmm.BIC;         % 2-component BIC
    gofTable.nSamples(gofRow) = numel(myData);   % needed for the count criterion
    % ------------------------------------------------------------------------
    gofTable.myData{gofRow}=myData;
    gofTable.mu1(gofRow)=mu(1);
    gofTable.mu2(gofRow)=mu(2);
    gofTable.sigma1(gofRow)=sigma(1);
    gofTable.sigma2(gofRow)=sigma(2);
    % Save final standard deviations.
    gofTable.std1(gofRow) = stdComponents(1);
    gofTable.std2(gofRow) = stdComponents(2);

    % Save targeted starting values.
    gofTable.startMu1(gofRow) = start.mu(1);
    gofTable.startMu2(gofRow) = start.mu(2);
    gofTable.startStd1(gofRow) = sqrt(start.Sigma(:,:,1));
    gofTable.startStd2(gofRow) = sqrt(start.Sigma(:,:,2));

    % Save whether the final fit retained the requested properties.
    gofTable.mu2AboveOne(gofRow) = finalMu2AboveOne;
    gofTable.std2Larger(gofRow) = finalStd2Larger;
    gofTable.weights1(gofRow)=weights(1);
    gofTable.weights2(gofRow)=weights(2);
    gofTable.ampUnweighted1(gofRow)=amp_unweighted(1);
    gofTable.ampUnweighted2(gofRow)=amp_unweighted(2);
    gofTable.ampWeighted1(gofRow)=amp_weighted(1);
    gofTable.ampWeighted2(gofRow)=amp_weighted(2);
    gofTable.boundaries{gofRow}=boundaries;
    gofTable.originalBounds{gofRow}=originalBounds;
    gofTable.minAllowedBoundary(gofRow)=boundaryThresh;
    if ~isempty(boundaries)
        gofTable.integration_log10(gofRow)=sum(myData(myData>boundaries));
        gofTable.integration(gofRow)=sum(LFPTable.LFP{gofRow}(LFPTable.LFP{gofRow}>10^boundaries));
        gofTable.percentIntegrated(gofRow)=sum(myData>boundaries)/length(myData)*100;
        gofTable.nMu1(gofRow)=sum(myData>mu(1)-2*sqrt(sigma(1)) & myData<mu(1)+2*sqrt(sigma(1)));
        gofTable.nMu2(gofRow)=sum(myData>mu(2)-2*sqrt(sigma(2)) & myData<mu(2)+2*sqrt(sigma(2)));
        gofTable.ttest_h(gofRow)=h;
        gofTable.ttest_p(gofRow)=p;
    else
        gofTable.integration_log10(gofRow)=0;
        gofTable.integration(gofRow)=0;
        gofTable.percentIntegrated(gofRow)=0;
        gofTable.nMu1(gofRow)=0;
        gofTable.nMu2(gofRow)=0;
    end
    
    % ===== EDIT 8: use figOutDir; sanitise filename (electrode names are fine,
    %               but region names contain '-' which is safe -- kept as-is).
    % saveas(gcf, fullfile("C:\BrewerLabResearch\Presentations&Figures\GaussFitFigsASR_HighGamma", "GMM " + axonNames(i) + ".png")) %Change b/w low and high gamma
    % ===== EDIT 8: use figOutDir.
    % OLD: saveas(gcf, fullfile(figOutDir, "GMM " + axonNames(i) + ".png"))
    saveas(gcf, fullfile(figOutDir, "GMM muSort " + axonNames(i) + ".png"))
    close(gcf)
    close(gcf)    
    % ======================================================================
    %gofRow=gofRow+1;
    % end
end

xlsxTable=gofTable;
xlsxTable.myData=[];
xlsxTable.gmm=[];
% ===== EDIT 9: originalBounds is a cell of VARIABLE-LENGTH vectors. writetable
%               rejects that as soon as any axon has 2 original crossings.
%               It worked so far only because every axon happened to have 1.
xlsxTable.originalBounds=[];
% ======================================================================
% OLD: xlsxTable.boundaries = cellfun(@(b) b(1), gofTable.boundaries);
% ===== NEW: b can now legitimately be [] -- cellfun(@(b)b(1)) errors on that.
bnd = nan(height(gofTable),1);
for k = 1:height(gofTable)
    if ~isempty(gofTable.boundaries{k}), bnd(k) = gofTable.boundaries{k}(1); end
end
xlsxTable.boundaries = bnd;
% ======================================================================
% writetable(xlsxTable,"C:\BrewerLabResearch\Presentations&Figures\GaussFitFigsASR_HighGamma\byAxonGMMFits.xlsx") %Change b/w low and high gamma
% OLD: writetable(xlsxTable, fullfile(figOutDir,"byAxonGMMFits.xlsx"))
writetable(xlsxTable, fullfile(figOutDir,"byAxonGMMFits_muSort.xlsx"))

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
%% Per-axon bimodality criteria

minWeight2    = 0.01;
minCount2     = 50;
minSeparation = 2;

% Only successfully fitted rows are eligible.
wasFitted = gofTable.fitStatus == "fitted";

% The two-component model must have the lower BIC.
bicPass = gofTable.BIC2 < gofTable.BIC1;

% The high-amplitude component must have sufficient fitted weight and
% sufficient estimated observations.
wPass = ...
    gofTable.weights2 >= minWeight2 & ...
    gofTable.weights2 .* gofTable.nSamples >= minCount2;

% The fitted means must be separated by at least minSeparation times the
% lower component's fitted standard deviation.
sepPass = ...
    gofTable.mu2 > ...
    gofTable.mu1 + minSeparation*sqrt(gofTable.sigma1);

% Confirm that a finite between-mean boundary exists.
hasBoundary = false(height(gofTable),1);

for r = 1:height(gofTable)
    b = gofTable.boundaries{r};

    hasBoundary(r) = ...
        ~isempty(b) && ...
        all(isfinite(b)) && ...
        isfinite(gofTable.mu2(r)) && ...
        gofTable.mu2(r) > b(1);
end

% Final non-destructive validity definition.
gofTable.mu2Valid = ...
    wasFitted & ...
    gofTable.mu2 > 1 & ...
    hasBoundary & ...
    bicPass & ...
    wPass & ...
    sepPass;

% Retain each individual criterion for inspection.
gofTable.bicPass = bicPass;
gofTable.wPass = wPass;
gofTable.sepPass = sepPass;
gofTable.hasBoundary = hasBoundary;

fprintf('\n--- mu2 validity breakdown (n = %d electrode pairs) ---\n', ...
    height(gofTable));

fprintf('  Successfully fitted : %d\n',sum(wasFitted));
fprintf('  Final mu2 > 1       : %d\n',sum(gofTable.mu2 > 1));
fprintf('  Has boundary        : %d\n',sum(hasBoundary));
fprintf('  BIC favors K=2      : %d\n',sum(bicPass));
fprintf('  Weight/count pass   : %d\n',sum(wPass));
fprintf('  Separation pass     : %d\n',sum(sepPass));
fprintf('  All criteria        : %d\n',sum(gofTable.mu2Valid));

for i = 1:length(all_reg)
    sel = gofTable.Subregion == all_reg(i) & hasLFP;

    fprintf('  %-10s %d/%d valid\n', ...
        all_reg(i), ...
        sum(sel & gofTable.mu2Valid), ...
        sum(sel));
end
% %% Per-axon bimodality criteria (REPLACES Dunnett testing and removal)
% %
% % --- WHY THE OLD CELL WAS REMOVED ---------------------------------------
% % The Dunnett test compared a control group of region mu1 values (n = #axons)
% % against one group per axon containing a SINGLE mu2 value (n = 1). With n=1
% % you cannot estimate within-group variance, so the pooled error term came
% % from the between-axon scatter in mu1. The rejection rule then deleted an
% % axon whenever p > 0.05, i.e. it used failure-to-reject from an underpowered
% % test as evidence of absence. That is why 4/5 regions were wiped out.
% %
% % Replaced with three per-axon criteria computed from each axon's own
% % thousands of amplitude samples. Non-destructive: writes a logical column
% % instead of nulling gofTable.boundaries, so re-running this cell is safe.
% % ------------------------------------------------------------------------
% 
% 
% minWeight2    = 0.01;   
% minCount2     = 50;     
% minSeparation = 2;      
% 
% bicPass = gofTable.BIC2 < gofTable.BIC1;                  % 2 comps beat 1 comp
% wPass   = gofTable.weights2 >= minWeight2 & ...
%           gofTable.weights2 .* gofTable.nSamples >= minCount2;
% sepPass = gofTable.mu2 > gofTable.mu1 + minSeparation*sqrt(gofTable.sigma1);
% 
% gofTable.mu2Valid = bicPass & wPass & sepPass;
% 
% gofTable.bicPass = bicPass;
% gofTable.wPass   = wPass;
% gofTable.sepPass = sepPass;
% 
% fprintf('\n--- mu2 validity breakdown (n = %d axons) ---\n', height(gofTable));
% fprintf('  BIC favours K=2 : %d\n', sum(bicPass));
% fprintf('  weight/count    : %d\n', sum(wPass));
% fprintf('  separation      : %d\n', sum(sepPass));
% fprintf('  ALL THREE       : %d\n', sum(gofTable.mu2Valid));
% for i = 1:length(all_reg)
%     sel = gofTable.Subregion == all_reg(i);
%     fprintf('  %-10s %d/%d valid\n', all_reg(i), ...
%             sum(sel & gofTable.mu2Valid), sum(sel));
% end
% 
% % ===== EDIT 11: per-region x per-criterion table. Read this to diagnose the
% %                missing orange bars before changing any threshold.
% fprintf('\n--- which criterion is pruning each region? ---\n');
% fprintf('%-10s %6s %6s %6s %6s %6s\n','region','n','BIC','wgt','sep','ALL');
% critTab = table();
% for i = 1:length(all_reg)
%     sel = gofTable.Subregion == all_reg(i);
%     fprintf('%-10s %6d %6d %6d %6d %6d\n', all_reg(i), sum(sel), ...
%             sum(sel & bicPass), sum(sel & wPass), sum(sel & sepPass), ...
%             sum(sel & gofTable.mu2Valid));
%     critTab.region(i)   = all_reg(i);
%     critTab.n(i)        = sum(sel);
%     critTab.bicPass(i)  = sum(sel & bicPass);
%     critTab.wPass(i)    = sum(sel & wPass);
%     critTab.sepPass(i)  = sum(sel & sepPass);
%     critTab.allPass(i)  = sum(sel & gofTable.mu2Valid);
% end
% writetable(critTab, fullfile(figOutDir,"mu2Validity_byRegion.xlsx"))
% 
% % Median separation actually achieved, in noise SDs -- tells you whether
% % minSeparation=2 is reasonable or arbitrarily harsh for this band.
% sepAchieved = (gofTable.mu2 - gofTable.mu1)./sqrt(gofTable.sigma1);
% fprintf('\nSeparation achieved (mu2-mu1)/SD1: median %.2f, 25th %.2f, 75th %.2f\n', ...
%         median(sepAchieved,'omitnan'), prctile(sepAchieved,25), prctile(sepAchieved,75));
% fprintf('Current minSeparation = %.2f\n', minSeparation);
% % ======================================================================
% 
% % ---- Per-region figures (same as before, now keyed off mu2Valid) --------
% ampRange=[min(cell2mat(LFPTable.LFP),[],"all"),max(cell2mat(LFPTable.LFP),[],"all")];
% binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
% x_grid = linspace(log10(ampRange(1)),log10(ampRange(2)),2000)';
% gauss = @(x, m, s) (1 / sqrt(2*pi*s)) .* exp(-((x - m).^2) / (2*s));
% 
% for i=1:length(all_reg)
%     figure
%     hold on
%     regRows=find(LFPTable.Subregion==all_reg(i));
%     goodMu2=[];
%     weight2=[];
% 
%     [maxAmp,~]=cellfun(@max,LFPTable.LFP(regRows));
%     [~,maxAmpAxonIdx]=max(maxAmp);
% 
%     for naxons=1:length(regRows)
%         r = regRows(naxons);
%         % --- OLD validity test (iscell branches + mu2 > boundaries) ---------
%         % if iscell(gofTable.boundaries(r)) ... mu2 > gofTable.boundaries{r} ...
%         % --- NEW ------------------------------------------------------------
%         if gofTable.mu2Valid(r)
%             goodMu2=[goodMu2,gofTable.mu2(r)];
%             weight2=[weight2,gofTable.weights2(r)];
%         end
%         % --------------------------------------------------------------------
% 
%         if naxons==maxAmpAxonIdx
%             histogram(log10(LFPTable.LFP{r}), log10(binEdges), 'Normalization','pdf', ...
%                       'FaceColor',[0.2 0.6 1.0], 'FaceAlpha',0.5)
%         end
%     end
% 
%     myMu1=gofTable.mu1(regRows);
%     mu1Avg = mean(myMu1);
%     newSD = std(myMu1);
%     meanPrior1=max(gofTable.weights1(regRows));
%     if isfinite(newSD) && newSD > 0
%         plot(x_grid,gauss(x_grid,mu1Avg,newSD).*meanPrior1,"LineWidth",4,"Color","b")
%     end
% 
%     if ~isempty(goodMu2)
%         mu2Avg = mean(goodMu2);
%         newSD = std(goodMu2);
%         if ~isfinite(newSD) || newSD == 0, newSD = eps; end
%         meanPrior2 = max(weight2);
%         plot(x_grid,gauss(x_grid,mu2Avg,newSD).*meanPrior2,"LineWidth",4,"Color","r")
%         area(x_grid,gauss(x_grid,mu2Avg,newSD).*meanPrior2, ...
%              'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none')
%     else
%         warning("%s: no valid mu2 to plot.", all_reg(i))
%     end
%     hold off
% 
%     ax=gca;
%     xticks([-2:4])
%     title(all_reg(i)+"  ("+sum(gofTable.mu2Valid(regRows))+"/"+numel(regRows)+" valid mu2)")
%     xlabel('Log_{10} \muV');
%     ylabel('Probability Density');
%     axis square
%     box on
%     set(ax,"FontSize",32)
%     ax.LineWidth=4;
%     ax.TickLength=[0.05,0.025];
%     xlim(log10(ampRange))
% 
%     saveas(gcf, fullfile(figOutDir, "RegionPDF " + all_reg(i) + ".png"))
%     % ==================================================================
% end
% 
% % ---- Export ------------------------------------------------------------
% xlsxTable=gofTable;
% xlsxTable.myData=[];
% xlsxTable.gmm=[];
% xlsxTable.originalBounds=[];                 % cell of variable-length vectors
% bnd = nan(height(gofTable),1);
% for k = 1:height(gofTable)
%     if ~isempty(gofTable.boundaries{k})
%         bnd(k) = gofTable.boundaries{k}(1);
%     end
% end
% xlsxTable.boundaries = bnd;
% % writetable(xlsxTable,"C:\BrewerLabResearch\Presentations&Figures\GaussFitFigsASR_HighGamma\byAxonGMMFits_post.xlsx") %Change b/w low and high gamma
% writetable(xlsxTable, fullfile(figOutDir,"byAxonGMMFits_post.xlsx"))   % ===== EDIT 12b

% %% ===== EDIT 13 (NEW CELL): Electrode-axon pair figures per subregion =====
% %
% %   PASS  -> thick solid line, legend entry marked (pass)
% %   FAIL  -> thin dotted line, legend entry marked (fail)
% %
% % Vertical lines:
% %   black dashed  = median GMM boundary across PASSING axons in the region
% %   grey dotted   = each passing axon's own GMM boundary (thin)
% %
% % Requires: gofTable.mu2Valid (from the cell above) and gofTable.boundaries.
% % Also uses regBinEdges/regFitVals saved in EDIT 3 if available.
% % =========================================================================
% 
% for i = 1:length(all_reg)
%     regRows = find(LFPTable.Subregion == all_reg(i));
%     if isempty(regRows)
%         warning('%s: no axons, skipping pair figure.', all_reg(i));
%         continue
%     end
% 
%     % per-region bin edges (reuse the histogram cell's if available)
%     if exist('regBinEdges','var') && numel(regBinEdges) >= i && ~isempty(regBinEdges{i})
%         bE = regBinEdges{i};
%     else
%         aR = [min(cell2mat(LFPTable.LFP(regRows)),[],"all"), ...
%               max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
%         aR(1) = max(aR(1), eps);
%         bE = logspace(log10(aR(1)), log10(aR(2)), 101);
%     end
%     centers = 10.^convert_edges_2_centers(log10(bE));
% 
%     pairColors = distinguishable_colors(numel(regRows), {'w','r','k'});
% 
%     figPair = figure('Color','w');
%     hold on
% 
%     legLabels = string.empty;
%     legHandles = [];
%     passBounds = [];
% 
%     for naxons = 1:numel(regRows)
%         r = regRows(naxons);
%         counts = histcounts(LFPTable.LFP{r}, bE);
% 
%         isPass = gofTable.mu2Valid(r);
%         if isPass
%             lsty = '-';  lwid = 4;
%             tag  = " (pass)";
%             thisBound = gofTable.boundaries{r};
%             if ~isempty(thisBound)
%                 passBounds = [passBounds, thisBound(1)];
%             end
%         else
%             lsty = ':';  lwid = 1.5;
%             tag  = " (fail)";
%         end
% 
%         hLine = plot(centers, counts, 'LineStyle', lsty, 'LineWidth', lwid, ...
%                      'Color', pairColors(naxons,:));
%         legHandles = [legHandles, hLine];
%         legLabels(end+1) = "FID " + LFPTable.fi(r) + " " + LFPTable.Electrode(r) + tag;
%     end
% 
%     % region-level Gauss fit (same red dashed + shaded area as your figure)
%     if exist('regFitVals','var') && numel(regFitVals) >= i && ~isempty(regFitVals{i})
%         fv = regFitVals{i};
%         hFit = plot(centers, fv, 'r--', 'LineWidth', 6);
%         area(centers, fv, 'FaceColor','r', 'FaceAlpha',0.25, 'EdgeColor','none');
%         legHandles = [legHandles, hFit];
%         legLabels(end+1) = "Gauss Fit";
%     end
% 
%     % each passing axon's own GMM boundary (thin grey), then the median (black)
%     for nb = 1:numel(passBounds)
%         xline(10^passBounds(nb), 'LineStyle',':', 'LineWidth',1, 'Color',[0.6 0.6 0.6]);
%     end
%     if ~isempty(passBounds)
%         medBound = median(passBounds);
%         xline(10^medBound, 'LineStyle','--', 'LineWidth',4, 'Color','k');
%     else
%         % no passing axons -> show the median of ALL boundaries for reference
%         allB = cellfun(@(b) (isempty(b))*NaN + (~isempty(b))*b(1)*(~isempty(b)), ...
%                        gofTable.boundaries(regRows));
%         if any(isfinite(allB))
%             xline(10^median(allB,'omitnan'), 'LineStyle','--', 'LineWidth',2, ...
%                   'Color',[0.5 0.5 0.5]);
%         end
%     end
% 
%     hold off
% 
%     axP = gca;
%     axP.XScale = "log";
%     xlabel("Amplitude (uV)")
%     ylabel("Count")
%     xticks(logspace(-3,4,8))
%     xticklabels("10^{"+[-3:4]+"}")
%     axis square
%     box on
%     axP.LineWidth = 2;
%     set(axP,"FontSize",14)
% 
%     lgd = legend(legHandles, legLabels, 'Location','eastoutside');
%     lgd.AutoUpdate = "off";
% 
%     nPass = sum(gofTable.mu2Valid(regRows));
%     text(0.025, 0.95, axP, all_reg(i) + " n=" + numel(regRows) + ...
%          "  (pass " + nPass + "/" + numel(regRows) + ")", ...
%          'Units','normalized', 'FontSize',14, 'FontWeight','bold');
% 
%     saveas(figPair, fullfile(figOutDir, "PairsGMMPass " + all_reg(i) + ".png"))
%     savefig(figPair, fullfile(figOutDir, "PairsGMMPass " + all_reg(i) + ".fig"))
% end
% 
% % ---- Companion summary figure: pass/fail count per region ---------------
% passCounts = zeros(length(all_reg),2);   % col1 = pass, col2 = fail
% for i = 1:length(all_reg)
%     sel = gofTable.Subregion == all_reg(i);
%     passCounts(i,1) = sum(sel & gofTable.mu2Valid);
%     passCounts(i,2) = sum(sel & ~gofTable.mu2Valid);
% end
% figSum = figure('Color','w');
% bPC = bar(reordercats(categorical(all_reg),all_reg), passCounts, 'stacked');
% legend({'pass GMM threshold','fail'},'Location','best')
% ylabel("Number of electrode-axon pairs")
% axis square; box on
% set(gca,"FontSize",14)
% saveas(figSum, fullfile(figOutDir, "PairPassFail_Summary.png"))

%% More Statistical tests and figure generation

if ~exist("gofTable","var")
   gofTable = readtable( ...
   fullfile(figOutDir,"byAxonGMMFits_muSort.xlsx"));
    % NOTE: this fallback cannot restore LFPTable (LFP is a cell column and
    % is not written to xlsx). You must re-run the Hilbert cell after a clear.
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
    validMu1Rows = regRows( ...
        gofTable.fitStatus(regRows) == "fitted" & ...
        isfinite(gofTable.mu1(regRows)));

    myMu1 = gofTable.mu1(validMu1Rows);

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
    validMu1Selection = ...
        sel & ...
        gofTable.fitStatus == "fitted" & ...
        isfinite(gofTable.mu1);

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