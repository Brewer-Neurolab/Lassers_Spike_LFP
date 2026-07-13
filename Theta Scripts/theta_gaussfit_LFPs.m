%% Fit and subtract gaussian per subregion

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

%% Stats calc
y=log10(lfpPropsTbl.ampIntegration);
g={lfpPropsTbl.Subregion};

[p,tbl,stats]=anovan(y,g);
[c,m]=multcompare(stats,'CriticalValueType','hsd');

figure
bar(regCats,m(:,1),'FaceColor','#75ebeb')
hold on
errorbar(regCats,m(:,1),m(:,2),'LineStyle','none','Color','k','LineWidth',3)
hold off
ylim([5,8])
% ylabel("Log_{10} Amplitude Integration")
ylabel("Mean Amplitude Integral x10^6")
yticks(0:1:10)
yticklabels(10.^(0:1:10)./10^6)
axis square
set(gca,"FontSize",28)
ax=gca;
ax.YAxis.MinorTick="on";
ax.YAxis.MinorTickValues=0:0.2:10;
box on
ax.LineWidth=4;
ax.TickLength=[0.05,0.025];

%% Anderson Darling testing

ADTable=LFPTable(:,1:6);
ADRow=1;
KurtosisTab=LFPTable(:,[1,2,4]);

for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    ampPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        pts=histogram(LFPTable.LFP{regRows(naxons)},binEdges);
        ampPts=pts.BinCounts;
        axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));

        set(gca,"XScale","log")

        LFPs=log10(LFPTable.LFP{regRows(naxons)}); %log transform data

        mu=fitVec{i}.b1;
        sigma=fitVec{i}.c1/sqrt(2);
        myGauss=makedist("Normal","mu",mu,"sigma",sigma);

        [h_gauss,p_gauss,adstat_gauss]=lillietest(LFPs);
        % [h_meanGauss,p_meanGauss,adstat_meanGauss]=lillietest(LFPs,'Distribution',myGauss);
        [h_EV,p_EV,adstat_EV]=lillietest(LFPs,'Distribution','extreme value');
        ADTable.h_gauss(ADRow)=h_gauss;
        ADTable.p_gauss(ADRow)=p_gauss;
        % ADTable.h_meanGauss(ADRow)=h_meanGauss;
        % ADTable.p_meanGauss(ADRow)=p_meanGauss;
        ADTable.h_EV(ADRow)=h_EV;
        ADTable.p_EV(ADRow)=p_EV;
        myKurtosis=kurtosis(LFPs);
        ADTable.Kurtosis(ADRow)=myKurtosis;

        KurtosisTab.Kurtosis(ADRow)=myKurtosis;

        ADRow=ADRow+1;

        figure
        plot(convert_edges_2_centers(log10(binEdges)),ampPts,'LineWidth',1)


        [f1,gof1]=fit(convert_edges_2_centers(log10(binEdges))',ampPts',"gauss1");
        hold on
        plot(f1,convert_edges_2_centers(log10(binEdges)),ampPts)
        hold off
        title(axonNames(naxons)+"Kurtosis="+myKurtosis)
        xlabel("Amplitude (\muV)")
        ylabel("Count")
    end
end

%% Single vs Double Gauss Fit Rsq BY AXON
close all
rng('default')
gofTable=LFPTable(:,1:6);
gofRow=1;

for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    binCenters=convert_edges_2_centers(log10(binEdges))';
    ampPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        hold on
        pts=histogram(log10(LFPTable.LFP{regRows(naxons)}),log10(binEdges));
        ampPts=pts.BinCounts;
        axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));

        [f1,gof1]=fit(convert_edges_2_centers(log10(binEdges))',ampPts',"gauss1");
        opts=fitoptions('Method','NonlinearLeastSquares');
        opts.StartPoint=[...
            10000, 0, 0.5,...
            2000,1.5,0.5];
        [f2,gof2]=fit(binCenters,ampPts',"gauss2",opts);
        [m1,v1]=lognstat(f1.b1,f1.c1);
        [m2,v2]=lognstat([f2.b1,f2.b2],[f2.c1,f2.c2]);

        % gofTable.Subregion(gofRow)=all_reg(i);
        gofTable.f1{gofRow}=f1;
        gofTable.gof1{gofRow}=gof1;
        gofTable.rsq1(gofRow)=gof1.rsquare;
        gofTable.mu(gofRow)=f1.b1;
        gofTable.sigma(gofRow)=f1.c1;
        gofTable.m(gofRow)=m1;
        gofTable.sd(gofRow)=sqrt(v1);
        gofTable.f2{gofRow}=f2;
        gofTable.gof2{gofRow}=gof2;
        gofTable.rsq2{gofRow}=gof2.rsquare;
        gofTable.mu2_1(gofRow)=f2.b1;
        gofTable.mu2_2(gofRow)=f2.b2;
        gofTable.sigma2_1(gofRow)=f2.c1;
        gofTable.sigma2_2(gofRow)=f2.c2;
        gofTable.m2_1(gofRow)=m2(1);
        gofTable.sd2_1(gofRow)=sqrt(v2(1));
        gofTable.m2_2(gofRow)=m2(2);
        gofTable.sd2_2(gofRow)=sqrt(v2(2));
        gofTable.gauss2better(gofRow)=gof2.rsquare>gof1.rsquare;
        gofTable.r2_diff(gofRow)=gof2.rsquare-gof1.rsquare;

        gofRow=gofRow+1;

        % plot(f1,convert_edges_2_centers(log10(binEdges)),ampPts)
        % plot(f2,convert_edges_2_centers(log10(binEdges)),ampPts)
        plot(binCenters,f1(binCenters),"LineWidth",4)
        plot(binCenters,f2(binCenters),"LineWidth",4)
        hold off
    end
end
%% Tail Index Test
close all
rng('default')
gofTable=LFPTable(:,1:6);
gofRow=1;
linTailIndex=[];
tailIndex=[];
for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    linBinEdges=linspace(ampRange(1),ampRange(2),101);
    binCenters=convert_edges_2_centers(log10(binEdges))';
    linBinCenters=convert_edges_2_centers(linBinEdges)';
    linAmpPts=[];
    ampPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        hold on
        pts=histogram(LFPTable.LFP{regRows(naxons)},linBinEdges);
        linAmpPts=pts.BinCounts;
        
        axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));
        
        %Find Tail Index
        % paramEsts=gpfit(LFPTable.LFP{regRows(naxons)});
        % kHat=paramEsts(1);
        % sigmaHat=paramEsts(2);
        % plot(linBinCenters,gppdf(linBinCenters,kHat,sigmaHat),"LineWidth",4)

        plot(linBinCenters,linAmpPts,"LineWidth",4)
        ax=gca;
        % ax.XScale="log";
        % ax.YScale="log";        
        hold off

        figure
        hold on
        pts=histogram(log10(LFPTable.LFP{regRows(naxons)}),log10(binEdges));
        ampPts=pts.BinCounts;
        plot(binCenters,ampPts,"LineWidth",4)
        hold off

        %Find Tail Index Past Median
        
        x=log10(LFPTable.LFP{regRows(naxons)});
        q=quantile(x,0.5);
        y=x(x>q)-q;
        paramEsts=gpfit(y);
        kHat=paramEsts(1);
        sigmaHat=paramEsts(2);
        figure
        hold on
        bins=linspace(0,log10(ampRange(2))-q,50);
        h=bar(bins,histc(y,bins)/(length(y)*.25),'histc');
        ygrid=linspace(0,1.1*max(y),100);
        plot(ygrid,gppdf(ygrid,kHat,sigmaHat))
        hold off
    end
end
%% Single vs Double Gauss Fit Rsq BY AXON SUBTRACTION
close all
rng('default')
gofTable=LFPTable(:,1:6);
gofRow=1;

for i=1:length(all_reg)

    regRows=find(gofTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    binCenters=convert_edges_2_centers(log10(binEdges))';
    ampPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        hold on
        pts=histogram(log10(LFPTable.LFP{regRows(naxons)}),log10(binEdges));
        ampPts=pts.BinCounts;
        axonNames(naxons)=gofTable.Subregion(regRows(naxons))+" FID "+gofTable.fi(regRows(naxons))+" "+gofTable.Electrode(regRows(naxons));

        %calculate Gaussian Mixed Model
        GMModel=fitgmdist(log10(LFPTable.LFP{regRows(naxons)}),2);
        figure

        % Fit Original Data
        opts=fitoptions('Method','NonlinearLeastSquares');
        opts.StartPoint=[10000,0,0.5];
        [f1,gof1,output1]=fit(binCenters,ampPts',"gauss1");
        plot(binCenters,f1(binCenters),"LineWidth",4)
        title(axonNames(naxons))
        hold off
        saveas(gcf,"C:\Users\lasss\Documents\Research\Brewer Lab work\All Presentations\Sam\Spike Field Coherence\Gauss Fit Figs\"+axonNames(naxons)+".png")
        
        %Subtract Fit 1 From distribution
        subAmpPts=ampPts'-f1(binCenters);
        subAmpPts(subAmpPts<0)=0;
        figure
        hold on
        % bar(binCenters,subAmpPts)
        h=histogram;
        h.BinEdges=log10(binEdges);
        h.BinCounts=subAmpPts';

        %Double fit Subtracted
        opts=fitoptions('Method','NonlinearLeastSquares');
        opts.StartPoint=[...
            3000, 0, 0.1,...
            2000,1.25,0.1];
        [f2,gof2,output2]=fit(binCenters,subAmpPts,"gauss2",opts);

        %Single fit subtracted
        opts=fitoptions('Method','NonlinearLeastSquares');
        opts.StartPoint=[2000,1.25,0.05];
        [f3,gof3,output3]=fit(binCenters,subAmpPts,"gauss1",opts);

        [m1,v1]=lognstat(f1.b1,f1.c1);
        [m2,v2]=lognstat([f2.b1,f2.b2],[f2.c1,f2.c2]);
        [m3,v3]=lognstat(f3.b1,f3.c1);

        % gofTable.Subregion(gofRow)=all_reg(i);
        gofTable.f1{gofRow}=f1;
        gofTable.gof1{gofRow}=gof1;
        gofTable.rsq1(gofRow)=gof1.rsquare;
        gofTable.residuals1{gofRow}=output1.residuals;
        gofTable.mu(gofRow)=f1.b1;
        gofTable.sigma(gofRow)=f1.c1;
        gofTable.m(gofRow)=m1;
        gofTable.sd(gofRow)=sqrt(v1);
        gofTable.fSub2{gofRow}=f2;
        gofTable.gofSub2{gofRow}=gof2;
        gofTable.rsqSub2{gofRow}=gof2.rsquare;
        gofTable.residuals2{gofRow}=output2.residuals;
        gofTable.muSub2_1(gofRow)=f2.b1;
        gofTable.muSub2_2(gofRow)=f2.b2;
        gofTable.sigmaSub2_1(gofRow)=f2.c1;
        gofTable.sigmaSub2_2(gofRow)=f2.c2;
        gofTable.mSub2_1(gofRow)=m2(1);
        gofTable.sdSub2_1(gofRow)=sqrt(v2(1));
        gofTable.mSub2_2(gofRow)=m2(2);
        gofTable.sdSub2_2(gofRow)=sqrt(v2(2));
        gofTable.fSub1{gofRow}=f3;
        gofTable.gofSub1{gofRow}=gof3;
        gofTable.rsqSub1{gofRow}=gof3.rsquare;
        gofTable.residuals3{gofRow}=output3.residuals;
        gofTable.muSub1_1(gofRow)=f3.b1;
        gofTable.sigmaSub1_1(gofRow)=f3.c1;
        gofTable.mSub1_1(gofRow)=m3(1);
        gofTable.sdSub1_1(gofRow)=sqrt(v3(1));

        %Integrate above mean+2SD for each, save values
        gofTable.thresh(gofRow)=f1.b1+(2*f1.c1);
        myVals=log10(LFPTable.LFP{regRows(naxons)});
        gofTable.valsAboveThresh{gofRow}=myVals(myVals>gofTable.thresh(gofRow));
        gofTable.integration(gofRow)=sum(gofTable.valsAboveThresh{gofRow});

        gofRow=gofRow+1;

        % plot(f1,convert_edges_2_centers(log10(binEdges)),ampPts)
        % plot(f2,convert_edges_2_centers(log10(binEdges)),ampPts)
        plot(binCenters,f2(binCenters),"LineWidth",4)
        plot(binCenters,f3(binCenters),"LineWidth",4)
        title(axonNames(naxons))
        xlabel("\muV")
        ylabel("Count")
        hold off
        saveas(gcf,"C:\Users\lasss\Documents\Research\Brewer Lab work\All Presentations\Sam\Spike Field Coherence\Gauss Fit Figs\"+"subtracted "+axonNames(naxons)+".png")
    end
end
%% Stats calc
y=log10(gofTable.integration);
g={gofTable.Subregion};

[p,tbl,stats]=anovan(y,g);
[c,m]=multcompare(stats,'CriticalValueType','hsd');

figure
bar(regCats,m(:,1),'FaceColor','#75ebeb')
hold on
errorbar(regCats,m(:,1),m(:,2),'LineStyle','none','Color','k','LineWidth',3)
hold off
ylim([5,8])
% ylabel("Log_{10} Amplitude Integration")
ylabel("Mean Amplitude Integral x10^6")
yticks(0:1:10)
yticklabels(10.^(0:1:10)./10^6)
axis square
set(gca,"FontSize",28)
ax=gca;
ax.YAxis.MinorTick="on";
ax.YAxis.MinorTickValues=0:0.2:10;
box on
ax.LineWidth=4;
ax.TickLength=[0.05,0.025];
%% Single vs Double Gauss Fit Rsq BY SUBREGION
close all
rng('default')
gofTable=table();
gofRow=1;

for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    binCenters=convert_edges_2_centers(log10(binEdges))';
    ampPts=[];
    axonNames=string();

    figure
    hold on
    pts=histogram(log10([LFPTable.LFP{regRows}]),log10(binEdges));
    ax=gca;
    % ax.XScale="log";
    ampPts=pts.BinCounts;
    % axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));

    [f1,gof1]=fit(binCenters,ampPts',"gauss1");

    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint=[...
        NaN, 0, NaN,...
        1*10^5,3,0.5];
    [f2,gof2]=fit(binCenters,ampPts',"gauss2",opts);

    [m1,v1]=lognstat(f1.b1,f1.c1);
    [m2,v2]=lognstat([f2.b1,f2.b2],[f2.c1,f2.c2]);

    gofTable.Subregion(gofRow)=all_reg(i);
    gofTable.f1{gofRow}=f1;
    gofTable.gof1{gofRow}=gof1;
    gofTable.rsq1(gofRow)=gof1.rsquare;
    gofTable.mu(gofRow)=f1.b1;
    gofTable.sigma(gofRow)=f1.c1;
    gofTable.m(gofRow)=m1;
    gofTable.sd(gofRow)=sqrt(v1);
    gofTable.f2{gofRow}=f2;
    gofTable.gof2{gofRow}=gof2;
    gofTable.rsq2{gofRow}=gof2.rsquare;
    gofTable.mu2_1(gofRow)=f2.b1;
    gofTable.mu2_2(gofRow)=f2.b2;
    gofTable.sigma2_1(gofRow)=f2.c1;
    gofTable.sigma2_2(gofRow)=f2.c2;
    gofTable.m2_1(gofRow)=m2(1);
    gofTable.sd2_1(gofRow)=sqrt(v2(1));
    gofTable.m2_2(gofRow)=m2(2);
    gofTable.sd2_2(gofRow)=sqrt(v2(2));
    gofTable.gauss2better(gofRow)=gof2.rsquare>gof1.rsquare;
    gofTable.r2_diff(gofRow)=gof2.rsquare-gof1.rsquare;

    gofRow=gofRow+1;

    % plot(f1,convert_edges_2_centers(log10(binEdges)),ampPts)
    % plot(f2,convert_edges_2_centers(log10(binEdges)),ampPts)

    plot(binCenters,f1(binCenters),"LineWidth",4)
    % plot(binCenters,f2(binCenters),"LineWidth",4)
    plot(binCenters,f2.a1.*exp(-((binCenters-f2.b1)./f2.c1).^2),"LineWidth",4)
    plot(binCenters,f2.a2.*exp(-((binCenters-f2.b2)./f2.c2).^2),"LineWidth",4)

    title(all_reg(i)+ " R^2 = "+round(gof1.rsquare,3)+", "+round(gof2.rsquare,3))
    subtitle(" ")
    hold off
    set(ax,"FontSize",40)
    axis square
    xlabel("Log_{10} Amplitude \muV")
    ylabel("Count")
    xlim([-2,4])
end

%% Sinusoid test
t=1/1000:1/1000:300;
% y=sin(t*10);
% y=pinknoise(length(t));
myGauss=makedist("Lognormal");
y=random(myGauss,[length(t),1]);

yAmps=abs(hilbert(y));

figure
plot(t,y)
hold on
plot(t,yAmps)
hold off

ampRange=[min(yAmps),max(yAmps)];
binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);

figure
hold on
pts=histogram(log10(yAmps),log10(binEdges));
ampPts=pts.BinCounts;

f=fit(convert_edges_2_centers(log10(binEdges))',ampPts',"gauss1");

plot(f,convert_edges_2_centers(log10(binEdges))',ampPts')
hold off

[h,p]=lillietest(log10(y));
figure
wblplot(log10(y))

%% OTSU THRESH
close all
otsuTable=LFPTable(:,1:6);
otsuRow=1;

for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    % axonsCounted=axonsCounted+length(regRows);
    % disp(axonsCounted+"/"+height(LFPTable))
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);
    ampPts=[];
    axonNames=string();
    for naxons=1:length(regRows)
        figure
        hold on
        pts=histogram(log10(LFPTable.LFP{regRows(naxons)}),log10(binEdges));
        ampPts=pts.BinCounts;
        axonNames(naxons)="FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));

        [T10,EM]=otsuthresh(ampPts);

        xline(T10)
        title(axonNames(naxons)+" T="+T10+" EM="+EM)
        xlabel("Amplitude (\muV)")
        ylabel("Count")
        hold off

        figure
        hold on
        pts=histogram(LFPTable.LFP{regRows(naxons)},binEdges);
        ampPts=pts.BinCounts;
        axonNames(naxons)=LFPTable.Subregion(regRows(naxons))+" FID "+LFPTable.fi(regRows(naxons))+" "+LFPTable.Electrode(regRows(naxons));

        xline(10^T10)
        title(axonNames(naxons)+" T="+10^T10+" EM="+EM)
        xlabel("Amplitude (\muV)")
        ylabel("Count")
        ax=gca;
        ax.XScale="log";
        hold off

        otsuTable.log10Thresh(otsuRow)=T10;
        otsuTable.log10Integration(otsuRow)=log10(sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=T10)));
        otsuTable.Thresh(otsuRow)=10^T10;
        otsuTable.Integration(otsuRow)=sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=(10^T10)));
        otsuTable.Effectiveness(otsuRow)=EM;
        otsuTable.PercentAboveThresh(otsuRow)=sum(LFPTable.LFP{regRows(naxons)}>=(10^T10))/length(LFPTable.LFP{regRows(naxons)})*100;

        otsuRow=otsuRow+1;
    end
end

% calculate log means
ampMeanTab=table();
for i=1:length(all_reg)
    ampMeanTab.subregion(i)=all_reg(i);
    ampMeanTab.mean(i)=mean(otsuTable.log10Integration(otsuTable.Subregion==all_reg(i) & otsuTable.Effectiveness>=0.65)); %& otsuTable.PercentAboveThresh>=0.20));
    ampMeanTab.sd(i)=std(otsuTable.log10Integration(otsuTable.Subregion==all_reg(i) & otsuTable.Effectiveness>=0.65));% & otsuTable.PercentAboveThresh>=0.20));
    ampMeanTab.min(i)=ampMeanTab.mean(i)-(2*ampMeanTab.sd(i));
    ampMeanTab.max(i)=ampMeanTab.mean(i)+(2*ampMeanTab.sd(i));
    ampMeanTab.n(i)=sum(otsuTable.Subregion==all_reg(i));
    % ampMeanTab.meanCount(i)=mean(otsuTable.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
    %remove outliers +/- 2SD
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & (lfpPropsTbl.ampIntegration<=10^ampMeanTab.min(i) | lfpPropsTbl.ampIntegration>=10^ampMeanTab.max(i)),:)=[];
    % remove integrations less than 50000
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.ampIntegration<=50000,:)=[];
    % remove time less than 5 seconds
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.sIntegration<=30,:)=[];
end

%
regCats=categorical(all_reg);
regCats=reordercats(regCats,all_reg);
figure
bar(regCats,ampMeanTab.mean)
hold on
errorbar(regCats,ampMeanTab.mean,ampMeanTab.sd./sqrt(ampMeanTab.n),'LineStyle','none','Color','k')
hold off
ylim([5,7])
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)

% calculate linear means
ampMeanTab=table();
for i=1:length(all_reg)
    ampMeanTab.subregion(i)=all_reg(i);
    ampMeanTab.mean(i)=mean(otsuTable.Integration(otsuTable.Subregion==all_reg(i) & otsuTable.Effectiveness>=0.65));% & otsuTable.PercentAboveThresh>=0.20));
    ampMeanTab.sd(i)=std(otsuTable.Integration(otsuTable.Subregion==all_reg(i) & otsuTable.Effectiveness>=0.65));% & otsuTable.PercentAboveThresh>=0.20));
    ampMeanTab.min(i)=ampMeanTab.mean(i)-(2*ampMeanTab.sd(i));
    ampMeanTab.max(i)=ampMeanTab.mean(i)+(2*ampMeanTab.sd(i));
    ampMeanTab.n(i)=sum(otsuTable.Subregion==all_reg(i));
    % ampMeanTab.meanCount(i)=mean(otsuTable.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
    %remove outliers +/- 2SD
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & (lfpPropsTbl.ampIntegration<=10^ampMeanTab.min(i) | lfpPropsTbl.ampIntegration>=10^ampMeanTab.max(i)),:)=[];
    % remove integrations less than 50000
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.ampIntegration<=50000,:)=[];
    % remove time less than 5 seconds
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.sIntegration<=30,:)=[];
end

%
regCats=categorical(all_reg);
regCats=reordercats(regCats,all_reg);
figure
bar(regCats,ampMeanTab.mean)
hold on
errorbar(regCats,ampMeanTab.mean,ampMeanTab.sd./sqrt(ampMeanTab.n),'LineStyle','none','Color','k')
hold off
% ylim([5,7])
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)

% Effectiveness Histogram
figure
histogram(otsuTable.Effectiveness,20)
axis square
xlabel("Otsu Effectiveness")
ylabel("Count")
ax=gca;
ax.FontSize=32;
xticks(0:0.05:1)

% Perc Above Thresh Histogram
figure
histogram(otsuTable.PercentAboveThresh,20)
axis square
xlabel("Percent Above Thresh")
ylabel("Count")
ax=gca;
ax.FontSize=32;
xticks(0:5:100)

%% Stats calc

% Test for Effectiveness Threshold
% y=log10(otsuTable.Integration(otsuTable.Effectiveness>=0.7));
% g={otsuTable.Subregion(otsuTable.Effectiveness>=0.7)};

%Test for time above thresh
% y=log10(otsuTable.Integration(otsuTable.PercentAboveThresh>=30));
% g={otsuTable.Subregion(otsuTable.PercentAboveThresh>=30)};

%Test for effectiveness and time
% y=log10(otsuTable.Integration(otsuTable.Effectiveness>=0.65 & otsuTable.PercentAboveThresh>=30));
% g={otsuTable.Subregion(otsuTable.Effectiveness>=0.65 & otsuTable.PercentAboveThresh>=30)};

%Test for thresh
y=log10(otsuTable.Integration(otsuTable.Thresh>=3.8));
g={otsuTable.Subregion(otsuTable.Thresh>=3.8)};

[p,tbl,stats]=anovan(y,g);
[c,m]=multcompare(stats,'CriticalValueType','hsd');

figure
bar(regCats,m(:,1),'FaceColor','#75ebeb')
hold on
errorbar(regCats,m(:,1),m(:,2),'LineStyle','none','Color','k','LineWidth',3)
hold off
ylim([4,8])
% ylabel("Log_{10} Amplitude Integration")
ylabel("Mean Amplitude Integral x10^6")
yticks(0:1:10)
yticklabels(10.^(0:1:10)./10^6)
axis square
set(gca,"FontSize",28)
ax=gca;
ax.YAxis.MinorTick="on";
ax.YAxis.MinorTickValues=0:0.2:10;
box on
ax.LineWidth=4;
ax.TickLength=[0.05,0.025];
%% OTSU THRESH Per Subregion
close all
otsuAllTable=table();
otsuRow=1;

for i=1:length(all_reg)

    regRows=find(LFPTable.Subregion==all_reg(i));
    ampRange=[min(cell2mat(LFPTable.LFP(regRows)),[],"all"),max(cell2mat(LFPTable.LFP(regRows)),[],"all")];
    binEdges=logspace(log10(ampRange(1)),log10(ampRange(2)),101);

    figure
    hold on
    regAmps=[LFPTable.LFP{regRows}];
    pts=histogram(log10(regAmps),log10(binEdges));
    ampPts=pts.BinCounts;

    [T10,EM]=otsuthresh(ampPts);

    xline(T10,'r--',"LineWidth",3)
    title(all_reg(i)+" T="+T10+" EM="+EM)
    subtitle(" ")
    xlabel("Log_{10} Amplitude (\muV)")
    ylabel("Count")
    ax=gca;
    ax.FontSize=32;
    axis square
    hold off

    figure
    hold on
    pts=histogram(regAmps,binEdges);
    ampPts=pts.BinCounts;

    xline(10^T10,'r--',"LineWidth",3)
    title(all_reg(i)+" T="+10^T10+" EM="+EM)
    subtitle(" ")
    xlabel("Amplitude (\muV)")
    ylabel("Count")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=32;
    axis square
    hold off

    otsuAllTable.subregion(otsuRow)=all_reg(i);
    otsuAllTable.log10Thresh(otsuRow)=T10;
    otsuAllTable.log10Integration(otsuRow)=log10(sum(regAmps(regAmps>=T10)));
    otsuAllTable.Thresh(otsuRow)=10^T10;
    otsuAllTable.Integration(otsuRow)=sum(regAmps(regAmps>=(10^T10)));
    otsuAllTable.Effectiveness(otsuRow)=EM;

    otsuRow=otsuRow+1;

end

%Individual Axon Integration
otsuIntegrationTab=LFPTable(:,1:6);
otsuRow=1;
for i=1:length(all_reg)
    regRows=find(LFPTable.Subregion==all_reg(i));
    for naxons=1:length(regRows)
        otsuIntegrationTab.log10Thresh(otsuRow)=otsuAllTable.log10Thresh(i);
        % otsuIntegrationTab.log10Integration(otsuRow)=log10(sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=otsuAllTable.log10Thresh(i))));
        otsuIntegrationTab.log10Integration(otsuRow)=log10(sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=log10(4))));
        otsuIntegrationTab.Thresh(otsuRow)=otsuAllTable.Thresh(i);
        % otsuIntegrationTab.Integration(otsuRow)=sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=(otsuAllTable.Thresh(i))));
        otsuIntegrationTab.Integration(otsuRow)=sum(LFPTable.LFP{regRows(naxons)}(LFPTable.LFP{regRows(naxons)}>=4));
        otsuRow=otsuRow+1;
    end
end

% calculate log means
ampMeanTab=table();
for i=1:length(all_reg)
    ampMeanTab.subregion(i)=all_reg(i);
    ampMeanTab.mean(i)=mean(otsuIntegrationTab.log10Integration(otsuIntegrationTab.Subregion==all_reg(i)));
    ampMeanTab.sd(i)=std(otsuIntegrationTab.log10Integration(otsuIntegrationTab.Subregion==all_reg(i)));
    ampMeanTab.min(i)=ampMeanTab.mean(i)-(2*ampMeanTab.sd(i));
    ampMeanTab.max(i)=ampMeanTab.mean(i)+(2*ampMeanTab.sd(i));
    ampMeanTab.n(i)=sum(otsuIntegrationTab.Subregion==all_reg(i));
    % ampMeanTab.meanCount(i)=mean(otsuTable.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
    %remove outliers +/- 2SD
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & (lfpPropsTbl.ampIntegration<=10^ampMeanTab.min(i) | lfpPropsTbl.ampIntegration>=10^ampMeanTab.max(i)),:)=[];
    % remove integrations less than 50000
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.ampIntegration<=50000,:)=[];
    % remove time less than 5 seconds
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.sIntegration<=30,:)=[];
end

%
regCats=categorical(all_reg);
regCats=reordercats(regCats,all_reg);
figure
bar(regCats,ampMeanTab.mean)
hold on
errorbar(regCats,ampMeanTab.mean,ampMeanTab.sd./sqrt(ampMeanTab.n),'LineStyle','none','Color','k')
hold off
ylim([5,7])
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)

% calculate linear means
ampMeanTab=table();
for i=1:length(all_reg)
    ampMeanTab.subregion(i)=all_reg(i);
    ampMeanTab.mean(i)=mean(otsuIntegrationTab.Integration(otsuIntegrationTab.Subregion==all_reg(i)));
    ampMeanTab.sd(i)=std(otsuIntegrationTab.Integration(otsuIntegrationTab.Subregion==all_reg(i)));
    ampMeanTab.min(i)=ampMeanTab.mean(i)-(2*ampMeanTab.sd(i));
    ampMeanTab.max(i)=ampMeanTab.mean(i)+(2*ampMeanTab.sd(i));
    ampMeanTab.n(i)=sum(otsuIntegrationTab.Subregion==all_reg(i));
    % ampMeanTab.meanCount(i)=mean(otsuTable.countIntegration(lfpPropsTbl.Subregion==all_reg(i)));
    %remove outliers +/- 2SD
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & (lfpPropsTbl.ampIntegration<=10^ampMeanTab.min(i) | lfpPropsTbl.ampIntegration>=10^ampMeanTab.max(i)),:)=[];
    % remove integrations less than 50000
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.ampIntegration<=50000,:)=[];
    % remove time less than 5 seconds
    % lfpPropsTbl(lfpPropsTbl.Subregion==all_reg(i) & lfpPropsTbl.sIntegration<=30,:)=[];
end

%
regCats=categorical(all_reg);
regCats=reordercats(regCats,all_reg);
figure
bar(regCats,ampMeanTab.mean)
hold on
errorbar(regCats,ampMeanTab.mean,ampMeanTab.sd./sqrt(ampMeanTab.n),'LineStyle','none','Color','k')
hold off
% ylim([5,7])
ylabel("Log_{10} Amplitude Integration")
axis square
set(gca,"FontSize",28)
%% Stats calc
% Test for Effectiveness Threshold
% y=log10(otsuIntegrationTab.Integration(otsuIntegrationTab.Effectiveness>=0.7));
% g={otsuIntegrationTab.Subregion(otsuIntegrationTab.Effectiveness>=0.7)};

%Test for time above thresh
% y=log10(otsuIntegrationTab.Integration(otsuIntegrationTab.PercentAboveThresh>=20));
% g={otsuIntegrationTab.Subregion(otsuIntegrationTab.PercentAboveThresh>=20)};

%Test for effectiveness and time
% y=log10(otsuIntegrationTab.Integration(otsuIntegrationTab.Effectiveness>=0.65 & otsuIntegrationTab.PercentAboveThresh>=30));
% g={otsuIntegrationTab.Subregion(otsuIntegrationTab.Effectiveness>=0.65 & otsuIntegrationTab.PercentAboveThresh>=30)};

y=log10(otsuIntegrationTab.Integration);
g={otsuIntegrationTab.Subregion};

[p,tbl,stats]=anovan(y,g);
[c,m]=multcompare(stats,'CriticalValueType','hsd');

figure
bar(regCats,m(:,1),'FaceColor','#75ebeb')
hold on
errorbar(regCats,m(:,1),m(:,2),'LineStyle','none','Color','k','LineWidth',3)
hold off
ylim([4,8])
% ylabel("Log_{10} Amplitude Integration")
ylabel("Mean Amplitude Integral x10^6")
yticks(0:1:10)
yticklabels(10.^(0:1:10)./10^6)
axis square
set(gca,"FontSize",28)
ax=gca;
ax.YAxis.MinorTick="on";
ax.YAxis.MinorTickValues=0:0.2:10;
box on
ax.LineWidth=4;
ax.TickLength=[0.05,0.025];
