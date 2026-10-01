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
    'Spikes','SpikeCount','Analytic_Sig','top_env','Phase',...
    'Power','AvgPower','BurstBounds','BurstBounds_avg_Power'}';

% initialize table
innit_tab=table('Size',[1,length(varnames)],...
    'VariableNames',varnames,...
    'VariableTypes',{'double','string','string','double','cell','cell',...
    'cell','cell','cell','cell','cell',...
    'cell','cell','cell','cell'});

%initialize structure to hold tables in
spike_amp_phase_struct.Delta=[];
spike_amp_phase_struct.Spindle=[];
spike_amp_phase_struct.Theta=[];
spike_amp_phase_struct.Low_Gamma=[];
spike_amp_phase_struct.High_Gamma=[];

%cut off frequencies
f1=[0.5,10,4,30,100]; %low cut off
f2=[4,16,10,100,300]; %high cut off

%bin sizes
freq_bins=1; %in seconds

% IMPORTANT NOTE: Allregion unit matched should have the same electrode
% pair order as the matching table, this should be automatic if using the
% correct preprocessing tools.
fs=25000;
re_fs=1000;

t=[0:1/fs:t_s-(1/fs)]; % 300 seconds
re_t=[0:1/re_fs:t_s-(1/re_fs)];

% Window Size For Frequency of Acceptable Samples
lfp_spike_win=[2,0.2,0.5,0.067,0.02]*re_fs; %seconds
lfp_median_thresh=[3,3,3,3,3];
lfp_median_step=0.5;

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

%load subregional spike trains
tunnel_parent_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
tunnel_folders=dir(tunnel_parent_dir);
tunnel_folder_names=string({tunnel_folders.name});
tunnel_folder_names=tunnel_folder_names([tunnel_folders.isdir]);
tunnel_folder_names=tunnel_folder_names(3:end);

% matching_table_tunnel=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\matching_table_wells_CW.mat");
% matching_table_tunnel=matching_table_tunnel.matching_table;

% Load synchrony measurements
% load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\coherence\AllSpikes_SpikeContrast_240530.mat")
% ff_axon_synchrony
%% downsample and filter wells
save_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\downsampled tunnels";
for fi=1:length(tunnel_folder_names)
    for subi=1:height(allregion_unit_matched_stim{fi})

        if length(allregion_unit_matched_stim{fi}.up_ff{subi})==1 & ...
            isempty(allregion_unit_matched_stim{fi}.up_fb{subi})
            elec2use=string(allregion_unit_matched_stim{fi}.("Electrode Pairs"){subi});
            elec2use=strsplit(elec2use,{'-'});
            elec2use=elec2use(1);
        elseif isempty(allregion_unit_matched_stim{fi}.up_ff{subi}) & ...
            length(allregion_unit_matched_stim{fi}.up_fb{subi})==1
            elec2use=string(allregion_unit_matched_stim{fi}.("Electrode Pairs"){subi});
            elec2use=strsplit(elec2use,{'-'});
            elec2use=elec2use(1);
        else
            continue
        end

        disp(elec2use+" "+string(fi))
        well_raw=load(fullfile(tunnel_parent_dir,tunnel_folder_names(fi),elec2use));
        well_raw=well_raw.data;
        re_LFP=resample(well_raw,re_fs,fs);
        sub_save_dir=fullfile(save_dir,"raw",tunnel_folder_names(fi));
        if ~isfolder(sub_save_dir)
            mkdir(sub_save_dir)
        end
        save(fullfile(sub_save_dir,elec2use+".mat"))
        for freqi=1:length(LFPs)
            % Filter out line noise
            lfreq=58;
            hfreq=62;
            fn=re_fs/2;
            Ws=[lfreq hfreq]/fn;
            [A,B,C,D]=butter(8,Ws,'stop');
            [sos,g] = ss2sos(A,B,C,D);
            re_LFP=filtfilt(sos,g,re_LFP);

            %silversmith filter
            [A,B,C,D] = butter(6,f1(freqi)/(re_fs/2),'high');
            [sos,g] = ss2sos(A,B,C,D);
            s_out1=filtfilt(sos,g,re_LFP);
            [A,B,C,D] = butter(8,f2(freqi)/(re_fs/2),'low');
            [sos,g] = ss2sos(A,B,C,D);
            filtered_data=filtfilt(sos,g,s_out1);

            sub_save_dir=fullfile(save_dir,LFPs{freqi},tunnel_folder_names(fi));
            if ~isfolder(sub_save_dir)
                mkdir(sub_save_dir)
            end
            save(fullfile(sub_save_dir,elec2use+".mat"))

        end
    end
    disp("Finished FID: "+string(fi))
end