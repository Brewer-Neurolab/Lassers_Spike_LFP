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
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

% all_reg=[interRegions,"EC-CA1"];
all_reg=[interRegions];

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
    nsamples_combine_thresh=(1/10)*re_fs*3;
    % nsamples_combine_thresh=[];

    %define min lfp length as 2x shortest theta cycle
    minLFPCycles=2; %default 2
    minLFPLength=(1/10)*minLFPCycles*re_fs;

    [LFPEndPts,LFPAmplitude,LFPHilbert]=identify_lfps(data,re_fs,t_rec,thresh_mult,minLFPLength,minLFPCycles,nsamples_combine_thresh);
    LFPAngles=wrapTo360(angle(LFPHilbert)*(180/pi));

    LFPTable.LFP{nFF}=LFPAmplitude;
    disp(nFF+"/"+height(ff_axon_tbl))
end

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
    regRows=find(LFPTable.Subregion==all_reg(i));
    axonsCounted=axonsCounted+length(regRows);
    disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    ampPts=[];
    cdfPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        pts=histogram(LFPTable.LFP{regRows(naxons)},binEdges);
        figure
        mycdf=histogram(LFPTable.LFP{regRows(naxons)},binEdges,"Normalization","cdf");
        ampPts(naxons,:)=pts.BinCounts;
        cdfPts(naxons,:)=mycdf.Values;
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

    last_bin=find(binEdges>m+2*v);
    last_bin=last_bin(1);

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

    regioncdf=histcounts([LFPTable.LFP{regRows}],binEdges,"Normalization","cdf");
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
ylim([5,7])
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)

%% GMM with Bayesian Decision Boundary — 1D
% -------------------------------------------------------
% INPUTS: Replace data_class1 and data_class2 with your
%         actual 1D data vectors (column vectors).
% -------------------------------------------------------

%% 1. Load / Define Your Data
close all
% Combine into one vector for GMM fitting
data_all = log10(LFPTable.LFP{regRows(naxons)})';

%% 2. Fit a 2-Component GMM Using EM
options = statset('MaxIter', 1000, 'TolFun', 1e-8);
gmm = fitgmdist(data_all, 2, 'Options', options, 'Replicates', 5);
% 'Replicates' runs EM from 5 different starting points
% to reduce the chance of landing on a local minimum

% Extract GMM Parameters
mu      = gmm.mu;                        % 2x1 means
sigma   = squeeze(gmm.Sigma);            % 2x1 variances (squeezed from 1x1x2)
weights = gmm.ComponentProportion;       % 1x2 mixture weights

fprintf('Component 1: mu = %.4f, sigma^2 = %.4f, weight = %.4f\n', ...
    mu(1), sigma(1), weights(1));
fprintf('Component 2: mu = %.4f, sigma^2 = %.4f, weight = %.4f\n', ...
    mu(2), sigma(2), weights(2));

%% 3. Define Class Priors
% Option A — Use mixture weights as priors (no prior knowledge)
prior1 = weights(1);
prior2 = weights(2);

% Option B — Use known class proportions if you have domain knowledge
% prior1 = 0.5; prior2 = 0.5;  % Uncomment to use equal priors

%% 4. Define Likelihoods and Posteriors
% Gaussian PDF for each component
gauss = @(x, m, s) (1 / sqrt(2*pi*s)) .* exp(-((x - m).^2) / (2*s));

% Evaluate over a fine grid for plotting and boundary finding
x_grid = linspace(min(data_all)-2, max(data_all)+2, 10000)';

% Likelihoods: P(x | class)
likelihood1 = gauss(x_grid, mu(1), sigma(1));
likelihood2 = gauss(x_grid, mu(2), sigma(2));

% Joint probabilities: P(x | class) * P(class)
joint1 = prior1 .* likelihood1;
joint2 = prior2 .* likelihood2;

% Evidence: P(x) — used to normalize posteriors
evidence = joint1 + joint2;

% Posteriors: P(class | x) via Bayes' theorem
posterior1 = joint1 ./ evidence;
posterior2 = joint2 ./ evidence;

%% 5. Find the Bayesian Decision Boundary
% The boundary is where posterior1 == posterior2, i.e., difference crosses 0
diff_posteriors = posterior1 - posterior2;
sign_changes = find(diff(sign(diff_posteriors)) ~= 0);

% There may be one or two crossing points depending on overlap
boundaries = x_grid(sign_changes);
fprintf('\nBayesian Decision Boundary/Boundaries:\n');
for i = 1:length(boundaries)
    fprintf('  Boundary %d: x = %.4f\n', i, boundaries(i));
end

%% 6. Classify Your Data Using the Boundary
% For a single boundary, classify based on side
if length(boundaries) == 1
    predicted_labels = double(data_all > boundaries(1)) + 1;
    % Class 1 = left of boundary, Class 2 = right of boundary
else
    % For two boundaries, assign based on nearest component mean
    predicted_labels = ones(size(data_all));
    for i = 1:length(data_all)
        p1 = prior1 * gauss(data_all(i), mu(1), sigma(1));
        p2 = prior2 * gauss(data_all(i), mu(2), sigma(2));
        if p2 > p1
            predicted_labels(i) = 2;
        end
    end
end

%% 7. Plot Results
figure;
hold on;

% Histogram of raw data
histogram(data_all, 40, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.6 1.0], 'FaceAlpha', 0.5, 'DisplayName', 'Class 1 Data');

% GMM component PDFs (scaled by prior)
plot(x_grid, joint1, 'b-', 'LineWidth', 2, 'DisplayName', 'GMM Component 1');
plot(x_grid, joint2, 'r-', 'LineWidth', 2, 'DisplayName', 'GMM Component 2');

% Overall GMM mixture PDF
plot(x_grid, evidence, 'k--', 'LineWidth', 1.5, 'DisplayName', 'GMM Mixture');

% Decision boundary line(s)
for i = 1:length(boundaries)
    xline(boundaries(i), 'm-', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Decision Boundary: x=%.3f', boundaries(i)));
end

xlabel('x');
ylabel('Probability Density');
title('1D GMM with Bayesian Decision Boundary');
legend('show', 'Location', 'best');
hold off;