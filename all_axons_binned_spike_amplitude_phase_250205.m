%% Spike vs Amplitude and Phase of LFPs
% Sam Lassers 11/20/23

%% Setup
clear all;
% close all;
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
load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat")

%load dynamics data
load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\spike_burst_dyn_table_stim.mat")
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
        load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\matching_table_cw.mat")
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
            cd("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels")
            cd(folders(j))

            temp_dir= convertStringsToChars(matching_table{i,2});
            disp("Processing FID: "+string(j)+" Tunnel: "+temp_dir)
            temp_dir(strfind(temp_dir,'-'):end)=[];
            Chanel_Name{i}={temp_dir};
            temp_dir1=[temp_dir,'.mat'];
            data=load (temp_dir1);
            re_LFP=data.re_LFP;
            %re_t=data.re_t;
            cd ..\..

            %load Spike from original
            cd("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A")
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

            for waves= 1%:length(LFPs)

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

%% by subregion plot of power distributuion
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\all\";

binEdges=logspace(-1,6,22);
binCenters=convert_edges_2_centers(binEdges);

pdf_table=table('Size',[0,5],'VariableTypes',{'double','string','string','string','cell'});
pdf_table.Properties.VariableNames=["FID","Subregion","Channel","FrequencyBand","PowerHistogram"];
% spike count v max amp

% for direction=[0,1] %0 is feed back, 1 is feedforward
for regi=1:length(subregions)
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    ax=gca;
    hold(ax, 'on')
    avg_power_count=[];
    se_power_count=[];
    for freq_bands=1%:length(LFPs)
        chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
        chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);
        my_chans=chan_power(:,[1,3]);
        chan_power=cell2mat(cellfun(@cell2mat, chan_power.Power, 'UniformOutput', false));


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
            h=histcounts(chan_power(chans,:),binEdges);
            f2=figure;
            hPDF=histogram(chan_power(chans,:),binEdges,'Normalization','probability');
            PDF_Vals=hPDF.Values;
            close gcf
            plot(ax2,binCenters,PDF_Vals,'LineWidth',3,'Color',my_colors(chans,:))
            % binned_power_mat=[binned_power_mat;h];
            binned_power_pdf=[binned_power_pdf;PDF_Vals];

            legend_log=[legend_log;"FID "+my_chans.FID(chans)+" "+my_chans.Channel(chans)];
        end
        lgd=legend(ax2,legend_log);
        xlabel("Power uV^2")
        ylabel("Power Probability")
        ax2.FontSize=30;
        xlim([binEdges(1),binEdges(end)])

        title(subregions(regi)+" All "+strrep(LFPs{freq_bands},'_',' '))
        % saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_pdf_avg.png")
        ax2.XScale="log";
        standardPlotParams_240422(ax2,"log","linear",[],[],0,[])
        set(ax2,"FontSize",24)
        saveas(gcf,save_dir+"\"+subregions(regi)+" "+LFPs{freq_bands}+"_ff_power_count_xlog_pdf_avg.png")
        % ax2.YScale="log";
        % saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_loglog_pdf_avg.png")

        hold(ax2,'off')
        avg_power_count=mean(binned_power_pdf,1,'omitmissing');
        se_power_count=std(binned_power_pdf,0,1,'omitmissing')/sqrt(size(binned_power_pdf,1));
        plot(ax,binCenters,avg_power_count,'LineWidth',3,'Color',band_colors(freq_bands,:))
        E=errorbar(ax,binCenters,avg_power_count,se_power_count,'LineStyle','none',...
            'Color','k','LineWidth',1);
        % set([E.bar,E.line],'truecoloralpha',255*0.5)
        pdf_table_temp=table();
        pdf_table_temp.FID=my_chans.FID;
        pdf_table_temp.Subregion=repmat(subregions(regi),[height(my_chans),1]);
        pdf_table_temp.Channel=my_chans.Channel;
        % pdf_table_temp.is_ff=repmat(direction,[height(my_chans),1]);
        pdf_table_temp.FrequencyBand=repmat(LFPs{freq_bands},[height(my_chans),1]);
        pdf_table_temp.PowerHistogram=num2cell(binned_power_pdf,2);
        pdf_table=[pdf_table;pdf_table_temp];
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
% end
%% K Means
kMeansTab=table();
row=1;
for regi=1:length(subregions)
    % f1=figure('units','normalized','outerposition',[0 0 1 1]);
    ax=gca;
    hold(ax, 'on')
    avg_power_count=[];
    se_power_count=[];
    for freq_bands=1%:length(LFPs)
        chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
        chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);
        my_chans=chan_power(:,[1,3]);
        chan_power=cell2mat(cellfun(@cell2mat, chan_power.Power, 'UniformOutput', false));
        % chan_power=reshape(chan_power,1,[]);
        kMeansTab.regi(row)=subregions(regi);
        kMeansTab.freq_band(row)=string(LFPs{freq_bands});
        kMeansTab.powers{row}=chan_power;
        [idx,C]=kmeans(chan_power,2);
        kMeansTab.idx{row}=idx;
        kMeansTab.centers{row}=C;
    end
end
%% Power integration over 10e2
close all
thresh=1*10^2;

powers=[];
power_regions=[];

for freq_bands=3%1:length(LFPs)
    ff_mean=[];
    fb_mean=[];
    ff_se=[];
    fb_se=[];
    ff_mean_logweight=[];
    fb_mean_logweight=[];
    ff_se_logweight=[];
    fb_se_logweight=[];
    ff_mean_weight=[];
    fb_mean_weight=[];
    ff_se_weight=[];
    fb_se_weight=[];
    for regi=1%:length(subregions)

        row2use=pdf_table.FrequencyBand==LFPs(freq_bands) & pdf_table.Subregion==subregions(regi);
        weights=cell2mat(pdf_table.PowerHistogram(row2use));
        weights=weights(:,binCenters>=thresh);
        pdf_table.threshSum(row2use)=sum(weights,2);
        pdf_table.threshSum_logweighted(row2use)=sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2);
        pdf_table.threshSum_weighted(row2use)=sum(weights.*repmat((binCenters(binCenters>=thresh)),[size(weights,1),1]),2);

        ff_mean=[ff_mean,mean(sum(weights,2))];
        ff_se=[ff_se,std(sum(weights,2))/sqrt(length(sum(weights,2)))];
        ff_mean_logweight=[ff_mean_logweight,mean(pdf_table.threshSum_logweighted(row2use))];
        ff_se_logweight=[ff_se_logweight,std(pdf_table.threshSum_logweighted(row2use))/sqrt(length(pdf_table.threshSum_logweighted(row2use)))];
        ff_mean_weight=[ff_mean_weight,mean(pdf_table.threshSum_weighted(row2use))];
        ff_se_weight=[ff_se_weight,std(pdf_table.threshSum_weighted(row2use))/sqrt(length(pdf_table.threshSum_weighted(row2use)))];

        powers=[powers;sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2);];
        power_regions=[power_regions;repmat(subregions(regi),size(sum(weights.*repmat(log10(binCenters(binCenters>=thresh)),[size(weights,1),1]),2)))];

    end
    myBars=[ff_mean_logweight];
    myError=[ff_se_logweight];
    myCat=categorical(subregions);
    myCat=reordercats(myCat,subregions);
    figure
    b=bar(myCat,myBars,'FaceColor','flat');
    hold on
    for k = 1:numel(b)                                                      % Recent MATLAB Versions
        xtips = b(k).XEndPoints;
        ytips = b(k).YEndPoints;
        errorbar(xtips,ytips,myError(:,k), '.k', 'MarkerSize',0.1)
    end

    hold off
    title(strrep(LFPs{freq_bands},"_"," "))
    ylabel("Average log signal power above "+string(thresh))
    set(gca,"FontSize",24)
end
%% Anova of powers above thresh
figure
[~,~,stats]=anovan(powers(powers>=0.25),{power_regions(powers>=0.25)});
[c,m]=multcompare(stats);
figure
bar(subregions,m(:,1))
hold on
errorbar(m(:,1),m(:,2),'LineStyle','none','Color','k')
hold off
%% by subregion plot of power distributuion PEAKS
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\all\";

binEdges=logspace(-1,4,16);
binCenters=convert_edges_2_centers(binEdges);

pdf_table=table('Size',[0,6],'VariableTypes',{'double','string','string','string','cell','cell'});
pdf_table.Properties.VariableNames=["FID","Subregion","Channel","FrequencyBand","PowerHistogram","ClusterCenters"];
% spike count v max amp

rng('default')
% for direction=[0,1] %0 is feed back, 1 is feedforward
for regi=1:length(subregions)
    f1=figure('units','normalized','outerposition',[0 0 1 1]);
    ax=gca;
    hold(ax, 'on')
    avg_power_count=[];
    se_power_count=[];
    for freq_bands=1%:length(LFPs)
        chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
        chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);
        my_chans=chan_power(:,[1,3]);
        chan_power=chan_power.pks;


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
            h=histcounts(chan_power{chans,:},binEdges);
            f2=figure;
            hPDF=histogram(chan_power{chans,:},binEdges,'Normalization','probability');
            PDF_Vals=hPDF.Values;
            close gcf
            plot(ax2,binCenters,PDF_Vals,'LineWidth',3,'Color',my_colors(chans,:))
            % binned_power_mat=[binned_power_mat;h];
            binned_power_pdf=[binned_power_pdf;PDF_Vals];

            legend_log=[legend_log;"FID "+my_chans.FID(chans)+" "+my_chans.Channel(chans)];
        end
        lgd=legend(ax2,legend_log);
        xlabel("Peak Amp uV")
        ylabel("Peak Amp Probability")
        ax2.FontSize=30;
        xlim([binEdges(1),binEdges(end)])

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

for freq_bands=1%:length(LFPs)
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
title("Delta")
ylabel("Log Amplitude (uV)")
set(gca,"FontSize",24)
ylim([1,3])
% axis square
%% by subregion plot of avg power distributuion
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\all\";

binEdges=logspace(-1,6,101);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp

for direction=[0,1] %0 is feed back, 1 is feedforward
    for regi=1:length(subregions)
        f1=figure('units','normalized','outerposition',[0 0 1 1]);
        ax=gca;
        hold(ax, 'on')
        avg_power_count=[];
        se_power_count=[];
        for freq_bands=1:length(LFPs)
            chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).is_ff==direction &...
                spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
            chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);

            chan_power=cell2mat(cellfun(@cell2mat, chan_power.Power, 'UniformOutput', false));

            % binned_power_mat=[];
            % binned_power_pdf=[];
            %plot all axons in one plot

            % hold on
            % for chans=1:size(chan_power,1)
            h=histcounts(chan_power,binEdges);
            f2=figure;
            hPDF=histogram(chan_power,binEdges,'Normalization','probability');
            PDF_Vals=hPDF.Values;
            close gcf
            % plot(binCenters,h,'LineWidth',3)
            % binned_power_mat=[binned_power_mat;h];
            binned_power_pdf=PDF_Vals;

            % end
            % avg_power_count=mean(binned_power_pdf,1,'omitmissing');
            se_bin_val=[];
            for nBin=1:length(binCenters)
                se_idx=find(chan_power>=binEdges(nBin) & chan_power<binEdges(nBin+1));
                se=std(chan_power(se_idx),0,"all")/sqrt(numel(se_idx));
                se_bin_val=[se_bin_val,se];
            end
            % se_power_count=std(chan_power,0,1,'omitmissing')/sqrt(numel(chan_power));
            plot(binCenters,binned_power_pdf,'LineWidth',3,'Color',band_colors(freq_bands,:))
            E=errorbar(binCenters,binned_power_pdf,se_bin_val,'LineStyle','none',...
                'Color','k','LineWidth',1);
            % set([E.bar,E.line],'truecoloralpha',255*0.5)
        end

        xlabel("Power uV^2")
        ylabel("Power Probability")
        ax.FontSize=30;
        xlim([binEdges(1),binEdges(end)])
        if direction==1
            title(subregions(regi)+" Feed-Forward")
            saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_pdf.png")
            ax.XScale="log";
            standardPlotParams_240422(ax,"log","linear",[],[],0,[])
            set(ax,"FontSize",24)
            saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_xlog_pdf.png")
            ax.YScale="log";
            saveas(gcf,save_dir+"\"+subregions(regi)+"_ff_power_count_loglog_pdf.png")
        else
            title(subregions(regi)+" Feedback")
            saveas(gcf,save_dir+"\"+subregions(regi)+"_fb_power_count_pdf.png")
            ax.XScale="log";
            standardPlotParams_240422(ax,"log","linear",[],[],0,[])
            set(ax,"FontSize",24)
            saveas(gcf,save_dir+"\"+subregions(regi)+"_fb_power_count_xlog_pdf.png")
            ax.YScale="log";
            saveas(gcf,save_dir+"\"+subregions(regi)+"_fb_power_count_loglog_pdf.png")
        end
        hold(ax, 'off')
    end
end
%% by frequency band plot of power distributuion
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\all\";

binEdges=logspace(-1,6,22);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp

for direction=[0,1] %0 is feed back, 1 is feedforward
    for freq_bands=1:length(LFPs)
        f1=figure('units','normalized','outerposition',[0 0 1 1]);
        ax=gca;
        hold(ax, 'on')
        avg_power_count=[];
        se_power_count=[];
        for regi=1:length(subregions)
            chans_log=spike_amp_phase_struct.(LFPs{freq_bands}).is_ff==direction &...
                spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi);
            chan_power=spike_amp_phase_struct.(LFPs{freq_bands})(chans_log,:);

            chan_power=cell2mat(cellfun(@cell2mat, chan_power.Power, 'UniformOutput', false));

            % binned_power_mat=[];
            binned_power_pdf=[];
            %plot all axons in one plot

            % hold on
            for chans=1:size(chan_power,1)
                h=histcounts(chan_power(chans,:),binEdges);
                f2=figure;
                hPDF=histogram(chan_power(chans,:),binEdges,'Normalization','probability');
                PDF_Vals=hPDF.Values;
                close gcf
                % plot(binCenters,h,'LineWidth',3)
                % binned_power_mat=[binned_power_mat;h];
                binned_power_pdf=[binned_power_pdf;PDF_Vals];

            end
            avg_power_count=mean(binned_power_pdf,1,'omitmissing');
            se_power_count=std(binned_power_pdf,0,1,'omitmissing')/sqrt(size(binned_power_pdf,1));
            plot(binCenters,avg_power_count,'LineWidth',3,'Color',band_colors(regi,:))
            E=errorbar(binCenters,avg_power_count,se_power_count,'LineStyle','none',...
                'Color','k','LineWidth',1);
            % set([E.bar,E.line],'truecoloralpha',255*0.5)
        end

        xlabel("Power uV^2")
        ylabel("Power Probability")
        ax.FontSize=30;
        xlim([binEdges(1),binEdges(end)])
        if direction==1
            title(strrep(LFPs{freq_bands},'_',' ')+" Feed-Forward")
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_ff_power_count_pdf.png")
            ax.XScale="log";
            standardPlotParams_240422(ax,"log","linear",[],[],0,[])
            set(ax,"FontSize",24)
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_ff_power_count_xlog_pdf.png")
            ax.YScale="log";
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_ff_power_count_loglog_pdf.png")
        else
            title(strrep(LFPs{freq_bands},'_',' ')+" Feedback")
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_fb_power_count_pdf.png")
            ax.XScale="log";
            standardPlotParams_240422(ax,"log","linear",[],[],0,[])
            set(ax,"FontSize",24)
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_fb_power_count_xlog_pdf.png")
            ax.YScale="log";
            saveas(gcf,save_dir+"\"+LFPs{freq_bands}+"_fb_power_count_loglog_pdf.png")
        end
        hold(ax, 'off')
    end
end

%% scatterplot spike count v max amplitude and spike count vs maxamp
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands}+"\scatter";
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        % figure
        % hold on
        % scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
        %     500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
        % xlabel("Max Amplitude (Micro Volts)")
        % ylabel("Spike Counts")
        % title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        % saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_scatter.png")
        % ax=gca;
        % ax.XScale="log";
        % ax.YScale="log";

        logMaxAmp=log(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels});
        % logMaxAmp=log(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels});
        logSpikeCount=log(spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels});

        %remove discontiuities
        okIdx=logSpikeCount>-Inf & logSpikeCount>1; %greater than 10 spikes for log10
        logMaxAmp=logMaxAmp(okIdx);
        logSpikeCount=logSpikeCount(okIdx);

        % const=polyfit(logMaxAmp,logSpikeCount,1);
        % m=const(1);
        % k=const(2);
        mdl=fitlm(logMaxAmp,logSpikeCount);
        k=mdl.Coefficients.Estimate(1);
        m=mdl.Coefficients.Estimate(2);
        % xfit=linspace(min(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),max(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),100);
        xfit=linspace(min(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}),max(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}),100);
        yfit=powerlawfitfun([k,m],xfit);
        % plot(xfit,yfit,"r--")
        % hold off
        figure
        plot(mdl)
        ylabel("log(Spike Count)")
        xlabel("log(Average Power uV^2)")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_scatter_loglog.png")
        hold off
        spike_amp_phase_struct.(LFPs{freq_bands}).FitSlopeAmp(channels)=m;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitInterceptAmp(channels)=k;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitRsqAmp(channels)=mdl.Rsquared.Ordinary;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitPValAmp(channels)=mdl.Coefficients.pValue(2);
    end
    %close all
end
%% scatterplot spike count v max amplitude and spike count vs AvgPower
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_AvgPower\tunnel\";

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands}+"\scatter";
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        % figure
        % hold on
        % scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
        %     500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
        % xlabel("Avg Power uV^2")
        % ylabel("Spike Counts")
        % title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        % saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
        %     spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_scatter.png")
        % ax=gca;
        % ax.XScale="log";
        % ax.YScale="log";

        % logMaxAmp=log(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels});
        logAvgAmp=log10(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels});
        logSpikeCount=log10(spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels});

        %remove discontiuities
        okIdx=logSpikeCount>-Inf & logSpikeCount>1; %greater than 10 spikes for log10
        logAvgAmp=logAvgAmp(okIdx);
        logSpikeCount=logSpikeCount(okIdx);

        % const=polyfit(logMaxAmp,logSpikeCount,1);
        % m=const(1);
        % k=const(2);
        mdl=fitlm(logAvgAmp,logSpikeCount);
        k=mdl.Coefficients.Estimate(1);
        m=mdl.Coefficients.Estimate(2);
        % xfit=linspace(min(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),max(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),100);
        xfit=linspace(min(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}),max(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}),100);
        yfit=powerlawfitfun([k,m],xfit);
        % plot(xfit,yfit,"r--")
        % hold off
        figure
        plot(mdl)
        ylabel("log(Spike Count)")
        xlabel("log(Average Power uV^2)")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_scatter_loglog.png")
        hold off
        spike_amp_phase_struct.(LFPs{freq_bands}).FitSlopePower(channels)=m;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitInterceptPower(channels)=k;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitRsqPower(channels)=mdl.Rsquared.Ordinary;
        spike_amp_phase_struct.(LFPs{freq_bands}).FitPValPower(channels)=mdl.Coefficients.pValue(2);
    end
    %close all
end
%% Compare pvalues and rsq for amp vs power
close all

rsqAvg=[];
pAvg=[];
rsqSE=[];
pSE=[];

for freq_bands=3
    for if_ff=[1,0]
        for regi=1:length(subregions)
            rows2use=spike_amp_phase_struct.(LFPs{freq_bands}).Subregion==subregions(regi) & spike_amp_phase_struct.(LFPs{freq_bands}).is_ff==if_ff;
            myRsqAmp=spike_amp_phase_struct.(LFPs{freq_bands}).FitRsqAmp(rows2use);
            myPAmp=spike_amp_phase_struct.(LFPs{freq_bands}).FitPValAmp(rows2use);
            myRsqPower=spike_amp_phase_struct.(LFPs{freq_bands}).FitRsqPower(rows2use);
            myPPower=spike_amp_phase_struct.(LFPs{freq_bands}).FitPValPower(rows2use);
            rsqAvg=[rsqAvg;[mean(myRsqAmp),mean(myRsqPower)]];
            pAvg=[pAvg;[mean(myPAmp),mean(myPPower)]];
            rsqSE=[rsqSE;[std(myRsqAmp)/sqrt(length(myRsqAmp)),std(myRsqPower)/sqrt(length(myRsqPower))]];
            pSE=[pSE;[std(myPAmp)/sqrt(length(myPAmp)),std(myPPower)/sqrt(length(myPPower))]];
        end
    end
end

x=categorical(subregions);
x=reordercats(x,subregions);
figure
b=bar(x,rsqAvg(1:5,:));
hold on
for k = 1:numel(b)                                                      % Recent MATLAB Versions
    xtips = b(k).XEndPoints;
    ytips = b(k).YEndPoints;
    errorbar(xtips,ytips,rsqSE(1:5,k), '.k', 'MarkerSize',0.1)
end
hold off
set(gca,"FontSize",18)
ylabel("Average Rsq")

x=categorical(subregions);
x=reordercats(x,subregions);
figure
b=bar(x,rsqAvg(6:10,:));
hold on
for k = 1:numel(b)                                                      % Recent MATLAB Versions
    xtips = b(k).XEndPoints;
    ytips = b(k).YEndPoints;
    errorbar(xtips,ytips,rsqSE(6:10,k), '.k', 'MarkerSize',0.1)
end
hold off
set(gca,"FontSize",18)
ylabel("Average Rsq")

x=categorical(subregions);
x=reordercats(x,subregions);
figure
b=bar(x,pAvg(1:5,:));
hold on
for k = 1:numel(b)                                                      % Recent MATLAB Versions
    xtips = b(k).XEndPoints;
    ytips = b(k).YEndPoints;
    errorbar(xtips,ytips,pSE(1:5,k), '.k', 'MarkerSize',0.1)
end
hold off
set(gca,"FontSize",18)
ylabel("Average P-Val")

x=categorical(subregions);
x=reordercats(x,subregions);
figure
b=bar(x,pAvg(6:10,:));
hold on
for k = 1:numel(b)                                                      % Recent MATLAB Versions
    xtips = b(k).XEndPoints;
    ytips = b(k).YEndPoints;
    errorbar(xtips,ytips,pSE(6:10,k), '.k', 'MarkerSize',0.1)
end
hold off
set(gca,"FontSize",18)
ylabel("Average P-Val")
%% plotting spike count v max amplitude histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands}+"\single channel histograms";
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        figure('units','normalized','outerposition',[0 0 1 1])
        %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
        %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
        spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
            spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

        hist_counts=[];
        for counts=2:length(binEdges)
            if counts~=length(binEdges)
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<binEdges(counts),1));
            else
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<=binEdges(counts),1));
            end
        end

        %h=bar(binCenters,hist_counts);
        plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Max Amplitude (Micro Volts)")
        ylabel("Spike Counts")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        ax.FontSize=30;

        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_maxamp.png")

        ax.XScale="log";
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_maxamp_xlog.png")
        ax.YScale="log";
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spikecount_v_maxamp_hist_loglog.png")

    end
    %close all
end
%% plotting spike count v max amplitude and spike count vs maxamp histograms zeros
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,101);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands}+"\zeros_histograms";
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        figure('units','normalized','outerposition',[0 0 1 1])
        %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
        %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
        spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
            spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

        hist_counts=[];
        for counts=2:length(binEdges)
            if counts~=length(binEdges)
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<binEdges(counts),1));
            else
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<=binEdges(counts),1));
            end
        end

        %h=bar(binCenters,hist_counts);
        plot(binCenters,hist_counts,'LineWidth',4)
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Max Amplitude (Micro Volts)")
        ylabel("Spike Counts")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        ax.FontSize=30;

        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_maxamp_zeros.png")

        ax.XScale="log";
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_maxamp_xlog_zeros.png")
        ax.YScale="log";
        saveas(gcf,save_dir_sub+"\log_"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spikecount_v_maxamp_hist_loglog_zeros.png")

    end
    %close all
end
%% one plot plotting spike count v max amplitude and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,101);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
            save_dir_sub=save_dir+"all"+LFPs{freq_bands}+"\histograms";
            if ~isfolder(save_dir_sub)
                mkdir(save_dir_sub)
            end

            if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                continue
            end

            %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
            %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
            %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
            spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

            hist_counts=[];
            for counts=2:length(binEdges)
                if counts~=length(binEdges)
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<binEdges(counts),1));
                else
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<=binEdges(counts),1));
                end
            end

            %h=bar(binCenters,hist_counts);
            plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4,...
                'Color',col_array_rgb(find(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)==subregions),:))
        end
        %close all
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Max Amplitude (Micro Volts)")
        ylabel("Spike Counts")
        ax.FontSize=30;

        if direction==1
            title(LFPs{freq_bands}+" Feed-Forward")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spikecount_v_maxamp_hist_loglog_all.png")
        else
            title(LFPs{freq_bands}+" Feedback")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spikecount_v_maxamp_hist_loglog_all.png")
        end

        hold off
    end
end
%% normalized one plot plotting spike count v max amplitude and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(-3,0,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
            save_dir_sub=save_dir+"all"+LFPs{freq_bands}+"\histograms";
            if ~isfolder(save_dir_sub)
                mkdir(save_dir_sub)
            end

            if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                continue
            end

            %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
            %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
            %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
            spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'./max(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels})];

            hist_counts=[];
            for counts=2:length(binEdges)
                if counts~=length(binEdges)
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<binEdges(counts),1));
                else
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<=binEdges(counts),1));
                end
            end

            %h=bar(binCenters,hist_counts);
            plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4,...
                'Color',col_array_rgb(find(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)==subregions),:))
        end
        %close all
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Max Amplitude (Micro Volts)")
        ylabel("Spike Counts")
        ax.FontSize=30;

        if direction==1
            title(LFPs{freq_bands}+" Feed-Forward")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_ff_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_ff_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_ff_spikecount_v_maxamp_hist_loglog_all.png")
        else
            title(LFPs{freq_bands}+" Feedback")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_fb_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_fb_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_normalized_fb_spikecount_v_maxamp_hist_loglog_all.png")
        end

        hold off
    end
end
%% one plot plotting spike count v max amplitude and spike count vs maxamp histograms zeros
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,101);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
            save_dir_sub=save_dir+"all"+LFPs{freq_bands}+"\zeros_histograms";
            if ~isfolder(save_dir_sub)
                mkdir(save_dir_sub)
            end

            if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                continue
            end

            %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
            %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
            %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
            spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

            hist_counts=[];
            for counts=2:length(binEdges)
                if counts~=length(binEdges)
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<binEdges(counts),1));
                else
                    hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                        & spike_v_maxamp(:,2)<=binEdges(counts),1));
                end
            end

            %h=bar(binCenters,hist_counts);
            plot(binCenters,hist_counts,'LineWidth',4,...
                'Color',col_array_rgb(find(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)==subregions),:))
        end
        %close all
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Max Amplitude (Micro Volts)")
        ylabel("Spike Counts")
        ax.FontSize=30;

        if direction==1
            title(LFPs{freq_bands}+" Feed-Forward")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_ff_spikecount_v_maxamp_hist_loglog_all.png")
        else
            title(LFPs{freq_bands}+" Feedback")
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spike_v_maxamp_all.png")
            ax.XScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spike_v_maxamp_xlog_all.png")
            ax.YScale="log";
            saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_fb_spikecount_v_maxamp_hist_loglog_all.png")
        end

        hold off
    end
end
%% by subregion plotting spike count v max amplitude and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        for regi=1:length(subregions)
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                save_dir_sub=save_dir+"all\"+LFPs{freq_bands}+"\histograms";
                if ~isfolder(save_dir_sub)
                    mkdir(save_dir_sub)
                end

                %skip if not right direction
                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                    continue
                end

                %skip if not right subregion
                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                    spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

                hist_counts=[];
                for counts=2:length(binEdges)
                    if counts~=length(binEdges)
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<binEdges(counts),1));
                    else
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<=binEdges(counts),1));
                    end
                end

                %h=bar(binCenters,hist_counts);
                plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
            end
            %close all
            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            xlabel("Max Amplitude (Micro Volts)")
            ylabel("Spike Counts")
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spikecount_v_maxamp_hist_loglog_all.png")
            else
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spikecount_v_maxamp_hist_loglog_all.png")
            end

            hold off
        end
    end
end
%% by subregion plotting spike count v max amplitude and spike count vs maxamp histograms zeros
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_maxamp\tunnel\";

binEdges=logspace(0,3.2,101);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        for regi=1:length(subregions)
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                save_dir_sub=save_dir+"all\"+LFPs{freq_bands}+"\zeros_histograms";
                if ~isfolder(save_dir_sub)
                    mkdir(save_dir_sub)
                end

                %skip if not right direction
                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                    continue
                end

                %skip if not right subregion
                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                    spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}'];

                hist_counts=[];
                for counts=2:length(binEdges)
                    if counts~=length(binEdges)
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<binEdges(counts),1));
                    else
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<=binEdges(counts),1));
                    end
                end

                %h=bar(binCenters,hist_counts);
                plot(binCenters,hist_counts,'LineWidth',4)
            end
            %close all
            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            xlabel("Max Amplitude (Micro Volts)")
            ylabel("Spike Counts")
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spikecount_v_maxamp_hist_loglog_all.png")
            else
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spikecount_v_maxamp_hist_loglog_all.png")
            end

            hold off
        end
    end
end

%% plotting spike count v avg POWER histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\";

binEdges=logspace(0,3.2,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands}+"\single channel histograms";
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        figure('units','normalized','outerposition',[0 0 1 1])
        %         scatter(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels},...
        %             500*ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels}),1),'.')
        spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
            spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}'];

        hist_counts=[];
        for counts=2:length(binEdges)
            if counts~=length(binEdges)
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<binEdges(counts),1));
            else
                hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                    & spike_v_maxamp(:,2)<=binEdges(counts),1));
            end
        end

        %h=bar(binCenters,hist_counts);
        plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
        xlim([min(binEdges),max(binEdges)])
        ax=gca;
        xlabel("Avgerage Power (uV^2)")
        ylabel("Spike Counts")

        if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)
            title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
                spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+" Feed-forward")
        else
            title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
                spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+" Feedback")
        end
        ax.FontSize=30;

        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_avgpower.png")

        ax.XScale="log";
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spike_v_avgpower_xlog.png")
        ax.YScale="log";
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_spikecount_v_avgpower_hist_loglog.png")

    end
    %close all
end

%% by subregion plotting spike count v average POWER and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\";

binEdges=logspace(-1,6,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        for regi=1:length(subregions)
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            nplots=1;
            p=[];
            mylegend="";
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                save_dir_sub=save_dir+"all\"+LFPs{freq_bands}+"\histograms";
                if ~isfolder(save_dir_sub)
                    mkdir(save_dir_sub)
                end

                %skip if not right direction
                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                    continue
                end

                %skip if not right subregion
                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                    spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}'];

                hist_counts=[];
                for counts=2:length(binEdges)
                    if counts~=length(binEdges)
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<binEdges(counts),1));
                    else
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<=binEdges(counts),1));
                    end
                end

                %h=bar(binCenters,hist_counts);
                %exclude zeros
                %plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
                %include zeros
                p(nplots)=plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4);
                mylegend(nplots)="FID " + string(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels)) + ...
                    " " + spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels);
                nplots=nplots+1;

            end
            %apply colors to plot
            mycolors=distinguishable_colors(nplots);
            legend(mylegend)

            %close all
            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            ax.ColorOrder=mycolors;
            xlabel("Mean Power (uV^2)")
            ylabel("Spike Counts")
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+strrep(LFPs{freq_bands},'_',' ')+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_spikecount_v_maxamp_hist_loglog_all.png")
            else
                title(subregions(regi)+" "+strrep(LFPs{freq_bands},'_',' ')+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_spikecount_v_maxamp_hist_loglog_all.png")
            end

            hold off
        end
    end
end
%% normalized by subregion plotting spike count v average POWER and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\";

binEdges=logspace(-3,0,41);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        for regi=1:length(subregions)
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                save_dir_sub=save_dir+"all\"+LFPs{freq_bands}+"\histograms";
                if ~isfolder(save_dir_sub)
                    mkdir(save_dir_sub)
                end

                %skip if not right direction
                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                    continue
                end

                %skip if not right subregion
                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                spike_v_maxamp=[spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels}',...
                    spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels}'./max(spike_amp_phase_struct.(LFPs{freq_bands}).AvgPower{channels})];

                hist_counts=[];
                for counts=2:length(binEdges)
                    if counts~=length(binEdges)
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<binEdges(counts),1));
                    else
                        hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                            & spike_v_maxamp(:,2)<=binEdges(counts),1));
                    end
                end

                %h=bar(binCenters,hist_counts);
                %exclude zeros
                %plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
                %include zeros
                plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
            end
            %close all
            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            xlabel("Mean Power (uV^2)")
            ylabel("Spike Counts")
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_ff_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_ff_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_ff_spikecount_v_maxamp_hist_loglog_all.png")
            else
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_fb_spike_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_fb_spike_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_normalized_fb_spikecount_v_maxamp_hist_loglog_all.png")
            end

            hold off
        end
    end
end
%% by subregion plotting BURST count v average POWER and spike count vs maxamp histograms
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_avgpower\tunnel\";

binEdges=logspace(-1,6,21);
binCenters=convert_edges_2_centers(binEdges);

% spike count v max amp
for freq_bands=3%1:length(LFPs)
    for direction=[0,1] %0 is feed back, 1 is feedforward
        hist_counts_all=[];
        for regi=1:length(subregions)
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                save_dir_sub=save_dir+"all\"+LFPs{freq_bands}+"\histograms";
                if ~isfolder(save_dir_sub)
                    mkdir(save_dir_sub)
                end

                %skip if not right direction
                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=direction
                    continue
                end

                %skip if not right subregion
                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                burst_v_maxamp=[ones(length(spike_amp_phase_struct.(LFPs{freq_bands}).BurstBounds_avg_Power{channels}'),1),...
                    spike_amp_phase_struct.(LFPs{freq_bands}).BurstBounds_avg_Power{channels}'];

                hist_counts=[];
                for counts=2:length(binEdges)
                    if counts~=length(binEdges)
                        hist_counts(counts-1)=sum(burst_v_maxamp(burst_v_maxamp(:,2)>=binEdges(counts-1)...
                            & burst_v_maxamp(:,2)<binEdges(counts),1));
                    else
                        hist_counts(counts-1)=sum(burst_v_maxamp(burst_v_maxamp(:,2)>=binEdges(counts-1)...
                            & burst_v_maxamp(:,2)<=binEdges(counts),1));
                    end
                end

                hist_counts_all=[hist_counts_all;hist_counts];
                %h=bar(binCenters,hist_counts);
                %exclude zeros
                %plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
                %include zeros
                plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
            end
            %close all
            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            xlabel("Mean Power (uV^2)")
            ylabel("Burst Counts")
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burst_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burst_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burstcount_v_maxamp_hist_loglog_all.png")
            else
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burst_v_maxamp_all.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burst_v_maxamp_xlog_all.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burstcount_v_maxamp_hist_loglog_all.png")
            end

            hold off

            % plot averages
            burst_avg=[];
            for nbins=1:length(binCenters)
                bin_counts=hist_counts_all(:,nbins);
                %bin_counts=bin_counts(bin_counts>0);
                bin_counts(bin_counts==0)=1;
                burst_avg(nbins)=mean(bin_counts);
                burst_std_err(nbins)=std(bin_counts)/sqrt(length(bin_counts));

                if isempty(bin_counts)
                    burst_avg(nbins)=0;
                    burst_std_err(nbins)=0;
                end
            end

            figure('units','normalized','outerposition',[0 0 0.5 1])
            hold on
            plot(binCenters(burst_avg>0),burst_avg(burst_avg>0),'Color','k','LineWidth',4)
            errorbar(binCenters(burst_avg>0),burst_avg(burst_avg>0),burst_std_err(burst_avg>0),...
                'LineStyle','none','LineWidth',2)
            hold off

            xlim([min(binEdges),max(binEdges)])
            ax=gca;
            xlabel("Mean Power (uV^2)")
            ylabel("Avg Burst/5 minutes")
            xlim([10^-1,10^6])
            xticks(logspace(-1,6,8))
            ax.FontSize=30;

            if direction==1
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feed-Forward")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burst_v_maxamp_avg.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burst_v_maxamp_xlog_avg.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_ff_burstcount_v_maxamp_hist_loglog_avg.png")
            else
                title(subregions(regi)+" "+LFPs{freq_bands}+" Feedback")
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burst_v_maxamp_avg.png")
                ax.XScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burst_v_maxamp_xlog_avg.png")
                ax.YScale="log";
                saveas(gcf,save_dir_sub+"\"+LFPs{freq_bands}+"_"+subregions(regi)+"_fb_burstcount_v_maxamp_hist_loglog_avg.png")
            end


        end
    end
end
%% plot spikes v phase -pi to pi
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_phase\tunnel\";
close all
%plot -pi to pi
phase_bin_edges=linspace(-pi,pi,21);
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        phasesig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Phase{channels});
        h=histogram(phasesig(spike_train),phase_bin_edges);
        xlabel("Phase")
        ylabel("Spike Count")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_phasehist_pi.png")
    end
end

%% plot spikes v phase -2pi to 2pi
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_phase\tunnel\";
close all
phase_bin_edges=linspace(-2*pi,2*pi,21*2);
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        phasesig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Phase{channels});

        %add 2pi and subtract 2 pi for signal
        phasesig_plus=phasesig+(2*pi);
        phasesig_minus=phasesig-(2*pi);

        h=histogram([phasesig(spike_train),phasesig_plus(spike_train),phasesig_minus(spike_train)],phase_bin_edges);
        xlabel("Phase")
        ylabel("Spike Count")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_phasehist_2pi.png")
    end
end

%% plot spikes v phase -360 to +360 degrees
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_phase\tunnel\";
close all
phase_bin_edges=linspace(-360,360,25);
for freq_bands=1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        phasesig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Phase{channels});

        %add 2pi and subtract 2 pi for signal
        phasesig=wrapTo360(phasesig*(180/pi));

        h=histogram([phasesig(spike_train),phasesig(spike_train)-360],phase_bin_edges);
        xlabel("Phase")
        ylabel("Spike Count")
        xticks([-360:30:360])
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))

        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_phasehist_360.png")
    end
end

%% Average Phase Bins
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_phase\tunnel\";
close all
phase_bin_edges=linspace(-360,360,25);
phase_bin_centers=convert_edges_2_centers(phase_bin_edges);
phase_ratios=[];
phase_ratios_pdf=[];
spike_counts_subreg=[];

for freq_bands=3%1:length(LFPs)
    for is_ff=[0,1]
        for regi=1:length(subregions)
            hist_counts=[];
            hist_pdf=[];
            figure('Position',[250,0,1200,1200])
            for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                %figure

                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=is_ff
                    continue
                end

                spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
                phasesig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Phase{channels});

                %add 2pi and subtract 2 pi for signal
                phasesig=wrapTo360(phasesig*(180/pi));

                h=histcounts([phasesig(spike_train),phasesig(spike_train)-360],phase_bin_edges);
                hist_counts=[hist_counts;h];

                h2=histogram([phasesig(spike_train),phasesig(spike_train)-360],phase_bin_edges,'Normalization','probability');
                hist_pdf=[hist_pdf;h2.Values];

                % title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
                %     spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))

            end
            save_dir_sub=save_dir+"all"+"\"+LFPs{freq_bands};
            if ~isfolder(save_dir_sub)
                mkdir(save_dir_sub)
            end
            hist_means=mean(hist_counts,1);
            hist_std=std(hist_counts,0,1);
            hist_SE=hist_std./sqrt(size(hist_counts,1));

            h_mean=histogram('BinEdges',phase_bin_edges,'BinCounts',hist_means);
            hold on
            errorbar(phase_bin_centers,hist_means,hist_SE,'LineStyle','none','Color','k')
            hold off

            xlabel("Phase")
            ylabel("Spike Count")
            xticks([-360:30:360])

            title(subregions(regi))

            ax=gca;
            ax.FontSize=30;

            saveas(gcf,save_dir_sub+"\"+subregions(regi)+"_"+is_ff+"_phasehist_360_avg.png")

            [max_bin,max_phase_loc]=max(hist_means);
            max_phase=phase_bin_centers(max_phase_loc);

            [min_bin,min_phase_loc]=min(hist_means);
            min_phase=phase_bin_centers(min_phase_loc);

            phase_ratios=[phase_ratios;{subregions(regi)},{is_ff},...
                {max_bin},{min_bin},{max_bin/min_bin},{max_phase+360},{min_phase+360}];

            params=[subregions(regi),is_ff,sum(hist_counts,'all')/2,max_phase]; %divide by 2 because the phase range is doubled

            %Probability
            hist_pdf(any(isnan(hist_pdf),2),:)=[];
            hist_means=mean(hist_pdf,1);
            hist_std=std(hist_pdf,0,1);
            hist_SE=hist_std./sqrt(size(hist_pdf,1));

            h_mean=histogram('BinEdges',phase_bin_edges,'BinCounts',hist_means);
            hold on
            errorbar(phase_bin_centers,hist_means,hist_SE,'LineStyle','none','Color','k')
            hold off

            xlabel("Phase")
            ylabel("Spike Probability")
            xticks([-360:30:360])

            title(subregions(regi))

            ax=gca;
            ax.FontSize=30;

            saveas(gcf,save_dir_sub+"\"+subregions(regi)+"_"+is_ff+"_phasehist_360_avg_PDF.png")

            [max_bin,max_phase_loc]=max(hist_means);
            max_phase=phase_bin_centers(max_phase_loc);

            [min_bin,min_phase_loc]=min(hist_means);
            min_phase=phase_bin_centers(min_phase_loc);

            phase_ratios_pdf=[phase_ratios_pdf;{subregions(regi)},{is_ff},...
                {max_bin},{min_bin},{max_bin/min_bin},{max_phase+360},{min_phase+360}];

            params=[params,max_phase,size(hist_counts,1)];
            spike_counts_subreg=[spike_counts_subreg;params];
        end
    end
end

%% Sawtooth
close all
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spikecount_v_phase\tunnel\";

%plot -pi to pi
phase_bin_edges=linspace(-pi,pi,21);
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        phasesig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Phase{channels});
        plot(t,phasesig)
        hold on
        plot(t(spike_train),phasesig(spike_train),'.','MarkerSize',18)
        xlabel("Time")
        ylabel("Phase")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_sawtooth.fig")
        hold off
    end
end

%% Signal bins and spikes overlayed on signal
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spike_signal_overlay\tunnel\";
%re_t=resample(t,1000,25000);
close all
for freq_bands=3%1:length(LFPs)
    for channels=1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure('units','normalized','outerposition',[0 0 1 1])
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        sig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).LFP{channels});
        %re_sig=resample(sig,1000,25000);
        plot(re_t,sig,'k')
        hold on
        plot(re_t(spike_train),sig(spike_train),'.')

        % plot bin splits
        window=re_fs*freq_bins; %in samples
        current_window=[1:window];
        for nbins=1:t_s/freq_bins
            xline(re_t(current_window(end)),'b')
            current_window=current_window+window;
        end

        xlabel("Time (s)")
        ylabel("Amplitude (uV)")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        ax=gca;
        ax.FontSize=24;
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_LFP_spike_overlay.fig")
        hold off
    end
end

%% Signal bins and spikes overlayed on POWER signal
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\spike_signal_overlay\tunnel\";
%re_t=downsample(t,25);
close all
for freq_bands=3%1:length(LFPs)
    for channels=23%1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
        figure('units','normalized','outerposition',[0 0 1 1])
        %         plot(spike_amp_phase_struct.(LFPs{freq_bands}).MaxAmp{channels},...
        %             spike_amp_phase_struct.(LFPs{freq_bands}).SpikeCount{channels})
        spike_train=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Spikes{channels});
        sig=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).Power{channels});

        %calculate and plot threshold
        [threshold,sig_rms]=rms_threshold_power(sig,2);

        %re_sig=resample(sig,1000,25000);
        plot(re_t,sig,'k')
        hold on
        yline(threshold)
        plot(re_t(spike_train),sig(spike_train),'.')

        % plot bin splits
        window=re_fs*freq_bins; %in samples
        current_window=[1:window];
        for nbins=1:t_s/freq_bins
            xline(re_t(current_window(end)),'b')
            current_window=current_window+window;
        end

        xlabel("Time (s)")
        ylabel("Power (µV^2)")
        title(spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+" "+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels))
        save_dir_sub=save_dir+folders(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels))+"\"+LFPs{freq_bands};
        if ~isfolder(save_dir_sub)
            mkdir(save_dir_sub)
        end
        ax=gca;
        ax.FontSize=24;
        saveas(gcf,save_dir_sub+"\"+spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).Channel(channels)+"_"+...
            spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)+"_power_spike_overlay.fig")
        hold off
    end
end

%% Cross correlation between binary spike train and LFP
max_lag=1000; %ms
max_lag_samples=(max_lag/1000)*re_fs;

%test on FID1 EC-DG F2
figure
%subplot(2,1,1)
hold on
test_sig=cell2mat(spike_amp_phase_struct.Theta.LFP{42});
test_spike_train=cell2mat(spike_amp_phase_struct.Theta.Spikes{42});
%plot(re_fs,test_sig,'k')
%plot(re_fs(test_spike_train),test_sig(test_spike_train),'.','MarkerSize',12)
hold off
[r,lags]=xcorr(test_sig,test_spike_train,max_lag_samples,'normalized');
%subplot(2,1,2)
stem(lags,r)
disp("got here!")
%%
fb_ec_dg=spike_amp_phase_struct.Theta(spike_amp_phase_struct.Theta.is_ff==0 & ...
    spike_amp_phase_struct.Theta.Subregion=="EC-DG",:);

%plot loop to test groups

figure('units','normalized','outerposition',[0 0 1 1])
ax=gca;
ax.XScale='log';
ax.YScale='log';
hold on
binEdges=logspace(-1,6,41);
binCenters=convert_edges_2_centers(binEdges);
for i=1:height(fb_ec_dg)
    spike_v_maxamp=[fb_ec_dg.SpikeCount{i}',...
        fb_ec_dg.AvgPower{i}'];

    hist_counts=[];
    for counts=2:length(binEdges)
        if counts~=length(binEdges)
            hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                & spike_v_maxamp(:,2)<binEdges(counts),1));
        else
            hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                & spike_v_maxamp(:,2)<=binEdges(counts),1));
        end
    end
    plot(binCenters(hist_counts>0),hist_counts(hist_counts>0),'LineWidth',4)
end
xlim([min(binEdges),max(binEdges)])
ax=gca;
xlabel("Mean Power (uV^2)")
ylabel("Spike Counts")
ax.FontSize=30;
hold off

%saveas(gcf,save_dir+"grouped single channel.png")

% with zeros to align
hist_struct=[];
binEdges=logspace(-1,6,41);
binCenters=convert_edges_2_centers(binEdges);
for i=1:height(fb_ec_dg)
    spike_v_maxamp=[fb_ec_dg.SpikeCount{i}',...
        fb_ec_dg.AvgPower{i}'];

    hist_counts=[];
    for counts=2:length(binEdges)
        if counts~=length(binEdges)
            hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                & spike_v_maxamp(:,2)<binEdges(counts),1));
        else
            hist_counts(counts-1)=sum(spike_v_maxamp(spike_v_maxamp(:,2)>=binEdges(counts-1)...
                & spike_v_maxamp(:,2)<=binEdges(counts),1));
        end
    end

    hist_struct.binCenters(i,:)=binCenters;
    hist_struct.histCounts(i,:)=hist_counts;
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=30;

%group 1
g1_avg=[];
g1_std=[];
g1_stderr=[];

for nbins=1:length(binCenters)
    col_vals=hist_struct.histCounts(:,nbins);
    col_vals=col_vals(col_vals>0);
    g1_avg(nbins)=mean(col_vals);
    g1_std(nbins)=std(col_vals);
    g1_stderr(nbins)=std(col_vals)./sqrt(length(col_vals));
end

plot(binCenters,g1_avg,'r','LineWidth',4)
errorbar(binCenters,g1_avg,g1_stderr)
xlabel("Power uV^2")
ylabel("Average Spike Count")

%% Axon LFP relation 2 target well spike

load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\optimization_table_tunnel2well.mat")
close all
xbins=logspace(-2,3,51);
for is_ff=[0,1] %0 is feed back, 1 is feedforward
    for regi=1:length(subregions)
        for freq_bands=1%1:length(LFPs)
            f1=figure('units','normalized','outerposition',[0 0 1 1]);
            ax=gca;
            hold(ax, 'on')
            for channels=16%1:height(spike_amp_phase_struct.(LFPs{freq_bands}))
                %figure

                if spike_amp_phase_struct.(LFPs{freq_bands}).Subregion(channels)~=subregions(regi)
                    continue
                end

                if spike_amp_phase_struct.(LFPs{freq_bands}).is_ff(channels)~=is_ff
                    continue
                end


                % Cycle eligible well channels
                target_subregion=strsplit(subregions(regi),{'-'});
                if is_ff
                    target_subregion=target_subregion(2);
                else
                    target_subregion=target_subregion(1);
                end
                myWells=matching_table_well.electrode(matching_table_well.subregion==target_subregion);

                hilbert_data=cell2mat(spike_amp_phase_struct.(LFPs{freq_bands}).top_env{channels});

                %for testing
                % threshold=3.5;

                myColors=distinguishable_colors(length(myWells));
                for subi=1:length(myWells)

                    %for testing
                    nLFP_Peaks=optimization_table.nPeaks(...
                        optimization_table.frequency_band=="Delta" &...
                        optimization_table.FID==3 &...
                        optimization_table.subregion=="DG-CA3" &...
                        optimization_table.tunnel=="L6" &...
                        optimization_table.is_ff==0 &...
                        optimization_table.threshold==3.5 &...
                        optimization_table.well==myWells(subi));

                    spike_train=load(fullfile(well_parent_dir,well_folder_names(spike_amp_phase_struct.(LFPs{freq_bands}).FID(channels)),...
                        myWells(subi)+"_spikes.mat"));
                    spike_idx=round((spike_train.index/1000)*fs);
                    spike_idx=round(remap(spike_idx,1,7500000,1,300000));
                    tunnel_amps_at_spikes=hilbert_data(spike_idx);
                    bin_vals=histcounts(tunnel_amps_at_spikes,xbins)/nLFP_Peaks;
                    plot(convert_edges_2_centers(xbins),bin_vals,'Color',myColors(subi,:))
                    ax.XScale='log';
                end
                legend(myWells)
            end
        end
        hold(ax, 'off')
        % close gcf
    end
end


%% Local Functions

% calculates rms and sets threshold to rms_multiplier times the rms of the
% input signal
function [threshold,sig_rms]=rms_threshold_power(powersig,rms_multiplier)
sig_rms=rms(powersig);
threshold=rms_multiplier*sig_rms;
end

function mapped_vals=remap(x, in_min, in_max, out_min, out_max)

mapped_vals=(x-in_min)*(out_max-out_min) / (in_max-in_min) + out_min;

end