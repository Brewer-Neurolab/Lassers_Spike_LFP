%% Setup
clear;
close all;
clc
% load matching_table_ccw.mat
folders=["4x 24571 210715 21div 210806_1_mat_files",...
    "4x 24572 210715 21div 210806_1_mat_files",...
    "4x 24573 210715 21div 210806_1_mat_files",...
    "4x 24574 210715 21div 210806_1_mat_files",...
    "4x 33152 210715 21div 210806_1_mat_files",...
    "4x 33168 210715 21div 210806_1_mat_files"];
s_folders=["4x 24571 210715 21div 210806_1_mat_files",...
    "4x 24572 210715 21div 210806_1_mat_files",...
    "4x 24573 210715 21div 210806_1_mat_files",...
    "4x 24574 210715 21div 210806_1_mat_files",...
    "4x 33152 210715 21div 210806_1_mat_files",...
    "4x 33168 210715 21div 210806_1_mat_files"];

t_s=300; %seconds, length of recording

varnames={'FID','Subregion','Channel','is_ff','LFP','MaxAmp',...
    'Spikes','SpikeCount','Analytic_Sig','top_env','bot_env','Phase',...
    'Power','AvgPower','BurstBounds','BurstBounds_avg_Power'}';

% initialize table
innit_tab=table('Size',[1,length(varnames)],...
    'VariableNames',varnames,...
    'VariableTypes',{'double','string','string','double','cell','cell',...
    'cell','cell','cell','cell','cell','cell',...
    'cell','cell','cell','cell'});

%initialize structure to hold tables in
spike_amp_phase_struct.Delta=[];
spike_amp_phase_struct.Spindle=[];
spike_amp_phase_struct.Theta=[];
spike_amp_phase_struct.Low_Gamma=[];
spike_amp_phase_struct.High_Gamma=[];

%cut off frequencies
f1=[0.5,10,4,30,100]; %low cut off
f2=[2,16,10,100,300]; %high cut off

%bin sizes
freq_bins=[1]; %in seconds

% IMPORTANT NOTE: Allregion unit matched should have the same electrode
% pair order as the matching table, this should be automatic if using the
% correct preprocessing tools.
fs=25000;
re_fs=1000;

t=[0:1/fs:t_s-(1/fs)]; % 300 seconds
re_t=[0:1/re_fs:t_s-(1/re_fs)];

%load matching table
load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A\allregion_unit_matched_cleaned.mat")

%load dynamics data
load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A\spike_burst_dyn_table_stim.mat")
%strip axon numbers from strings
axon_names=spike_burst_dyn_table_stim.channel_name;
for i=1:length(axon_names)
    axon=axon_names(i);
    axon=strsplit(axon,{'-'});
    axon_names(i)=axon(1);
end
spike_burst_dyn_table_stim.channel_name=axon_names;

LFPs=fieldnames(spike_amp_phase_struct);

Chanel_Name=[];
is_cw=[1,1,1,1,1,1];
qcw=[];

col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E","#7d0505"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; [125, 5, 5]./255];

subregions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];

powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1));

band_colors=[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
%% Filter for freq band and plot

%plot upstream waves

tic
for j=1:length(allregion_unit_matched_stim)
    if is_cw(j)==1
        load("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A\4x 24571 210715 21div 210806_1_mat_files\matching_table_wNMI.mat")
        qcw=[qcw,'T'];
    else
        load matching_table_ccw.mat
        qcw=[qcw,'F'];
    end
    Interval(1:height(matching_table))=cell(1,height(matching_table));
    Length(1:height(matching_table))=cell(1,height(matching_table));
    Integral_conv(1:height(matching_table))=cell(1,height(matching_table));

    row=1;

    for i=1:21
        %select for tunnels with LFP
        if ~isempty(allregion_unit_matched_stim{j}.ff_cdt{i}) | ~isempty(allregion_unit_matched_stim{j}.fb_cdt{i})

            %load file from down samples
            cd("C:\BrewerLabResearch\OneDrive_1_7-16-2025\downsampled tunnels\High_Gamma")
            cd(folders(j))

            temp_dir= convertStringsToChars(matching_table{i,2});
            disp("Processing FID: "+string(j)+" Tunnel: "+temp_dir)
            temp_dir(strfind(temp_dir,'-'):end)=[];
            Chanel_Name{i}={temp_dir};
            temp_dir1=[temp_dir,'.mat'];
            data=load (temp_dir1);
            re_LFP=data.filtered_data;
            %re_t=data.re_t;
            cd ..\..

            %load Spike from original
            cd("C:\BrewerLabResearch\OneDrive_1_7-16-2025\18-Apr-2023_A")
            cd(s_folders(j))
            temp_dir= convertStringsToChars(matching_table{i,2});
            temp_dir(strfind(temp_dir,'-'):end)=[];
            Chanel_Name{i}={temp_dir};
            temp_dir2=strcat('times_',temp_dir,'.mat');
            temp_dir=[temp_dir,'.mat'];
            data=load (temp_dir);
            %data2=data.data;

            try
                data=load (temp_dir2, 'cluster_class');
                spike_data=data.cluster_class  (:,2);
                spike_data=spike_data(data.cluster_class(:,1)==1);
            catch
                spike_data=[];
            end
            cd ..\..

            % Select only a portion of the data, modulate using t_s in s
            re_LFP=re_LFP(1:(t_s*re_fs));
            %data2=data2(1:(t_s*fs));
            spike_data=spike_data/1000; %ms to s
            spike_data=spike_data(spike_data<=t_s);
            re_t=re_t(1:(t_s*re_fs));

            %round spike data for 1000 hZ
            spike_data_rounded=round(spike_data,3);

            for waves= 5%:length(LFPs)

                LFP_tab=innit_tab;

                % Filter out line noise
                lfreq=58;
                hfreq=62;
                fn=re_fs/2;
                Ws=[lfreq hfreq]/fn;
                [A,B,C,D]=butter(8,Ws,'stop');
                [sos,g] = ss2sos(A,B,C,D);
                re_LFP=filtfilt(sos,g,re_LFP);

                %filter
                %[A,B,C,D] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass'); % 3
                [A,B,C,D] = butter(6,f1(waves)/(re_fs/2),'high');
                %[B,A] = butter(3,[fc1/(fs/2),fc2/(fs/2)],'bandpass');
                [sos,g] = ss2sos(A,B,C,D);
                s_out1=filtfilt(sos,g,re_LFP);
                [A,B,C,D] = butter(8,f2(waves)/(re_fs/2),'low');
                [sos,g] = ss2sos(A,B,C,D);
                filtered_data=filtfilt(sos,g,s_out1);

                window=re_fs*freq_bins; %in samples

                current_window=[1:window];

                analitic_wins=[];
                maxamp_wins=[];
                phase_wins=[];
                spike_wins=[];
                spike_counts=[];
                LFP_win=[];
                power_win=[];
                avgpower_win=[];
                top_env_win=[];
                bot_env_win=[];

                %[logical_spikes]=ismembertol(t,spike_data,1e-10);
                [logical_spikes]=ismembertol(re_t,spike_data_rounded,1e-10);

                hilbert_data=hilbert(filtered_data);
                [top_env,bot_env]=envelope(filtered_data);

                for nbins=1:t_s/freq_bins
                    binned_spikes=logical_spikes(current_window);
                    sigphase = ((angle(hilbert_data(current_window))));
                    analiticsig = imag(hilbert_data(current_window));

                    %assign to arrays
                    %changed max amp to hilbert max amp 240430
                    %changed power to hilbert data squared
                    LFP_win=[LFP_win,{filtered_data(current_window)}];
                    spike_wins=[spike_wins,{binned_spikes}];
                    spike_counts=[spike_counts,sum(binned_spikes)];
                    analitic_wins=[analitic_wins,{analiticsig}];
                    maxamp_wins=[maxamp_wins,max(abs(hilbert_data(current_window)))];
                    phase_wins=[phase_wins,{sigphase}];
                    power_win=[power_win,{abs(hilbert_data((current_window))).^2}];
                    avgpower_win=[avgpower_win,mean(abs(hilbert_data((current_window))).^2)];
                    top_env_win=[top_env_win,{top_env(current_window)}];
                    bot_env_win=[bot_env_win,{bot_env(current_window)}];

                    current_window=current_window+window;

                end
                LFP_tab.FID(1)=j;
                LFP_tab.Subregion(1)=allregion_unit_matched_stim{j}.Subregion{i};
                LFP_tab(1,:).Channel=Chanel_Name{i};
                if length(allregion_unit_matched_stim{j}.ff_cdt{i})==1 && isempty(allregion_unit_matched_stim{j}.fb_cdt{i})
                    LFP_tab.is_ff(1)=1;
                else
                    LFP_tab.is_ff(1)=0;
                end

                LFP_tab.LFP{1}=LFP_win;
                LFP_tab.Spikes{1}=spike_wins;
                LFP_tab.SpikeCount{1}=spike_counts;
                LFP_tab.Analytic_Sig{1}=analitic_wins;
                LFP_tab.MaxAmp{1}=maxamp_wins;
                LFP_tab.Phase{1}=phase_wins;
                LFP_tab.Power{1}=power_win;
                LFP_tab.AvgPower{1}=avgpower_win;
                LFP_tab.top_env{1}=top_env_win;
                LFP_tab.bot_env{1}=bot_env_win;

                % Finds peaks of hilbert transform and adds height and position
                % to table within 1 cycle of relevant freq band

                [pks,locs]=findpeaks(top_env,'MinPeakHeight',2*mean(top_env),'MinPeakDistance',1/f2(waves)*re_fs);
                LFP_tab.pks{1}=pks;
                LFP_tab.pk_locs{1}=locs;

                % burstBounds=spike_burst_dyn_table_stim.BurstBounds{...
                %     spike_burst_dyn_table_stim.fi==j & spike_burst_dyn_table_stim.channel_name==Chanel_Name{i}};
                % %resample burstBounds
                % burstBounds=remap(burstBounds,min(1),max(length(t)),min(1),max(length(re_t)));
                % burstBounds=round(burstBounds);
                % LFP_tab.BurstBounds{1}=burstBounds;

                %get average power of signal within burst
                % for bursts=1:length(burstBounds)
                %     burst_idx=ismembertol(re_t,re_t(burstBounds(bursts,1):burstBounds(bursts,2)));
                %     avg_power_burst(bursts)=mean(filtered_data(burst_idx).^2);
                % end
                %
                % LFP_tab.BurstBounds_avg_Power{1}=avg_power_burst;

                % concatonates table to structure
                if isempty(spike_amp_phase_struct.(LFPs{waves}))
                    spike_amp_phase_struct.(LFPs{waves})=LFP_tab;
                else
                    spike_amp_phase_struct.(LFPs{waves})=[spike_amp_phase_struct.(LFPs{waves});LFP_tab];
                end

                row=row+1;
            end
            % remove top row because this table construction is a
            % little janky
        end
    end


    disp("Files complete: "+string(j))
end
toc

%save("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spike_amp_phase_struct.mat",'spike_amp_phase_struct')
%% reload
% too large to do this
% load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spike_amp_phase_struct.mat")
%% by subregion plot of power distributuion PEAKS
close all
save_dir="C:\BrewerLabResearch\OneDrive_1_7-16-2025\coherence\spikecount_v_avgpower\tunnel\all\High_Gamma";

binEdges=logspace(-1,4,16);
binCenters=convert_edges_2_centers(binEdges);

pdf_table=table('Size',[0,6],'VariableTypes',{'double','string','string','string','cell','cell'});
pdf_table.Properties.VariableNames=["FID","Subregion","Channel","FrequencyBand","PowerHistogram","ClusterCenters"];
% spike count v max amp

rng('default')

allintegrals = [];
regionCounts = [];
avgintegrals = zeros(1,length(subregions));
% for direction=[0,1] %0 is feed back, 1 is feedforward
for regi=1:length(subregions)
    f1=figure('units','normalized','outerposition',[0 0 1 1]);
    ax=gca;
    hold(ax, 'on')
    avg_power_count=[];
    se_power_count=[];
    numChans = 0;
    %%For each region of graph, find integral by doing:
    %%abs(bc2-bc1)*abs(Pval2-Pval1)/2 + abs(bc2-bc1) * min(Pval1, Pval2)
    for freq_bands=5%:length(LFPs)
        chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
        chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);
        my_chans=chan_power(:,[1,3]);
        chan_power=chan_power.pks;
        subintegrals = zeros(1,length(chan_power));

        % binned_power_mat=[];
        binned_power_pdf=[];
        %plot all axons in one plot

        % hold on
        figure('units','normalized','outerposition',[0 0 1 1]);
        ax2=gca;
        hold(ax2,'on')
        legend_log=[];
        my_colors=distinguishable_colors(size(chan_power,1));
        for chans=1:size(chan_power,1)
            currentInt = 0;
            h=histcounts(chan_power{chans,:},binEdges);
            f2=figure;
            hPDF=histogram(chan_power{chans,:},binEdges,'Normalization','probability');
            PDF_Vals=hPDF.Values;
            for li=1:length(PDF_Vals)
                if (li ~= 1 && PDF_Vals(li) ~= PDF_Vals(li-1))
                    currentInt=currentInt+(abs(binCenters(li+1)-binCenters(li))*(abs(PDF_Vals(li+1)-PDF_Vals(li))/2)) + (abs(binCenters(li+1)-binCenters(li)) * min(PDF_Vals(li), PDF_Vals(li+1)));
                end
            end
            subintegrals(chans) = currentInt;
            numChans = numChans + 1;
            close gcf
            plot(ax2,binCenters,PDF_Vals,'LineWidth',3,'Color',my_colors(chans,:))
            % binned_power_mat=[binned_power_mat;h];
            binned_power_pdf=[binned_power_pdf;PDF_Vals];

            legend_log=[legend_log;"FID "+my_chans.FID(chans)+" "+my_chans.Channel(chans)];
        end
        allintegrals = [allintegrals, subintegrals];
        regionCounts = [regionCounts, numChans];
        avgintegrals(regi) = mean(subintegrals);
        lgd=legend(ax2,legend_log);
        xlabel("Peak Amp uV")
        ylabel("Peak Amp Probability")
        ax2.FontSize=30;
        xlim([1,1000])

        title(subregions(regi)+" All "+strrep(LFPs{freq_bands},'_',' '))
        % saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_pdf_avg.png")
        ax2.XScale="log";
        standardPlotParams_240422(ax2,"log","linear",[],[],0,[])
        set(ax2,"FontSize",24)
        fontsize(lgd,10,"points")
        saveas(gcf,save_dir+"\"+subregions(regi)+" "+LFPs{freq_bands}+"_ff_power_count_xlog_pdf_avg.png")
        % ax2.YScale="log";
        % saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_loglog_pdf_avg.png")

        hold(ax2,'off')
        avg_power_count=mean(binned_power_pdf,1,'omitmissing');
        se_power_count=std(binned_power_pdf,0,1,'omitmissing')/sqrt(size(binned_power_pdf,1));
        plot(ax,binCenters,avg_power_count,'LineWidth',3,'Color',band_colors(freq_bands,:))
        E=errorbar(ax,binCenters,avg_power_count,se_power_count,'LineStyle','none',...
            'Color','k','LineWidth',1);
        axis square
        % set([E.bar,E.line],'truecoloralpha',255*0.5)

        % cluster to determine integration cutoff

        [idx,C,sumd,D]=kmeans(binned_power_pdf*log10(binCenters'),2);

        figure('units','normalized','outerposition',[0 0 1 1]);
        p3=plot(binCenters,binned_power_pdf');
        [p3(idx==find(C==max(C))).Color]=deal([1,0,0]);
        [p3(idx==find(C==min(C))).Color]=deal([0,0,1]);
        ax3=gca;
        ax3.XScale='log';
        lgd=legend(ax3,legend_log);
        title(subregions(regi)+" All "+strrep(LFPs{freq_bands},'_',' '))
        xlabel("Peak Amp uV")
        ylabel("Peak Amp Probability")
        ax3.FontSize=30;
        fontsize(lgd,10,"points")

        xline(10.^mean(C),"LineWidth",5,"LineStyle","--","Color","k")
        lgd.String=lgd.String(1:end-1);

        pdf_table_temp=table();
        pdf_table_temp.FID=my_chans.FID;
        pdf_table_temp.Subregion=repmat(subregions(regi),[height(my_chans),1]);
        pdf_table_temp.Channel=my_chans.Channel;
        % pdf_table_temp.is_ff=repmat(direction,[height(my_chans),1]);
        pdf_table_temp.FrequencyBand=repmat(LFPs{freq_bands},[height(my_chans),1]);
        pdf_table_temp.PowerHistogram=num2cell(binned_power_pdf,2);
        pdf_table_temp.ClusterCenters=repmat({C},height(my_chans),1);
        pdf_table=[pdf_table;pdf_table_temp];

        % add xline to colored plot
        xline(ax2,10.^mean(C),"LineWidth",5,"LineStyle","--","Color","k")
    end

    xlabel(ax,"Power uV^2")
    ylabel(ax,"Power Probability")
    ax.FontSize=30;
    xlim(ax,[binEdges(1),binEdges(end)])
    title(ax,subregions(regi))
    fontsize(lgd,10,"points")
    saveas(ax,save_dir+"\"+subregions(regi)+"_power_count_pdf_avg.png")
    ax.XScale="log";
    standardPlotParams_240422(ax,"log","linear",[],[],0,[])
    set(ax,"FontSize",24)
    saveas(ax,save_dir+"\"+subregions(regi)+"_power_count_xlog_pdf_avg.png")
    ax.YScale="log";
    saveas(ax,save_dir+"\"+subregions(regi)+"_power_count_loglog_pdf_avg.png")

    hold(ax, 'off')

end
pdf_tab_reload=pdf_table;
%% Power integration difference of centroid threshold
% close all
% thresh=1*10^2;
pdf_table=pdf_tab_reload;
myThresh=[];

powers=[];
power_regions=[];

for freq_bands=5%:length(LFPs)
    ff_mean=[];
    ff_se=[];
    ff_mean_logweight=[];
    ff_se_logweight=[];
    ff_mean_weight=[];
    ff_se_weight=[];
    for regi=1:length(subregions)

        row2use=pdf_table.FrequencyBand==LFPs(freq_bands) & pdf_table.Subregion==subregions(regi);
        thresh=pdf_table.ClusterCenters(row2use);
        thresh=10.^mean(thresh{1});
        myThresh=[myThresh,thresh];
        weights=cell2mat(pdf_table.PowerHistogram(row2use));
        weights=weights(:,binCenters>=thresh);
        pdf_table.threshSum(row2use)=sum(weights,2);
        pdf_table.threshSum_logweighted(row2use)=sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2);
        pdf_table.threshSum_weighted(row2use)=sum(weights.*repmat((binCenters(binCenters>=thresh)),[size(weights,1),1]),2);

        % pdf_table(pdf_table.threshSum==0,:)=[];
        % pdf_table(pdf_table.threshSum_logweighted==0,:)=[];
        % pdf_table(pdf_table.threshSum_weighted==0,:)=[];
        %
        % ff_mean=[ff_mean,mean(sum(weights,2))];
        % ff_se=[ff_se,std(sum(weights,2))/sqrt(length(sum(weights,2)))];
        % ff_mean_logweight=[ff_mean_logweight,mean(pdf_table.threshSum_logweighted(row2use))];
        % ff_se_logweight=[ff_se_logweight,std(pdf_table.threshSum_logweighted(row2use))/sqrt(length(pdf_table.threshSum_logweighted(row2use)))];
        % ff_mean_weight=[ff_mean_weight,mean(pdf_table.threshSum_weighted(row2use))];
        % ff_se_weight=[ff_se_weight,std(pdf_table.threshSum_weighted(row2use))/sqrt(length(pdf_table.threshSum_weighted(row2use)))];
        %
        % powers=[powers;sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2);];
        % power_regions=[power_regions;repmat(subregions(regi),size(sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2)))];

    end
    % myBars=[ff_mean_logweight];
    % myError=[ff_se_logweight];
    % myCat=categorical(subregions);
    % myCat=reordercats(myCat,subregions);
    % figure
    % b=bar(myCat,myBars,'FaceColor','flat');
    % hold on
    % for k = 1:numel(b)                                                      % Recent MATLAB Versions
    %     xtips = b(k).XEndPoints;
    %     ytips = b(k).YEndPoints;
    %     errorbar(xtips,ytips,myError(:,k), '.k', 'MarkerSize',0.1)
    % end
    %
    % hold off
    % title(strrep(LFPs{freq_bands},"_"," "))
    % ylabel("Average log signal power above threshold")
    % set(gca,"FontSize",24)

end
pdf_table(pdf_table.threshSum==0,:)=[];
%% Anova of powers above thresh
figure
[~,~,stats]=anovan(pdf_table.threshSum_logweighted./pdf_table.threshSum,{[pdf_table.Subregion+pdf_table.FrequencyBand]});
[c,m]=multcompare(stats);
figure
bar(subregions,m(:,1))
hold on
errorbar(m(:,1),m(:,2),'LineStyle','none','Color','k')
hold off

% TODO need to put in loop later
title("Average Amplitude (uV) Per Subregion")
xlabel("Subregion")
ylabel("Log Amplitude (uV)")
set(gca,"FontSize",24)
ylim([1,3])
% axis square

%% Histogram of log amplitude vs average amplitude of each subregion

logAmps = log10(avgintegrals);
allTracks = 1;
hold on
passNums = [];
for hi = 1:5
    numPass = 0;
    for ai = 1:regionCounts(hi)
        if (allintegrals(allTracks) >= myThresh(hi))
            numPass = numPass + 1;
        end
        allTracks = allTracks + 1;
    end
    passNums = [passNums, numPass];
    % %{text(x, y, labelText, ...
    %     'HorizontalAlignment', 'center', ...
    %     'VerticalAlignment', 'middle', ...
    %     'Color', 'w', ...
    %     'FontWeight', 'bold', ...
    %     'Rotation', 90); 
    % 
end

hold off

histogram('Categories',logAmps,'BinCounts',5)
axis square;
