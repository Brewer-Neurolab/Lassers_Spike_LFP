%% Xcorr CA1 Subregion to G10 spikes FID 5
clear
clc
G10=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\times_G12.mat");
G10=G10.cluster_class(:,2);
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

CA1_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==4 & well_spike_dyn.fi==5);

fs=25000;
t=1/fs:1/fs:300;
%% xcorr
xcorr_table=table();

G10_spike_train=zeros(1,length(t));
G10_spike_train(ismembertol(t,G10/1000))=1;

for nElec=1:length(CA1_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\"+CA1_Elec(nElec)+"_spikes.mat");
    my_spike_train=zeros(1,length(t));
    my_spike_train(ismembertol(t,myElec.index/1000))=1;
    [r,l]=xcorr(G10_spike_train,my_spike_train,fs*0.01,"normalized");
    figure
    plot(l/fs*1000,r)
    title(CA1_Elec(nElec))

    [~,P]=corrcoef(G10_spike_train,my_spike_train);

    xcorr_table.elec(nElec)=CA1_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.l{nElec}=l;
    xcorr_table.pval(nElec)=P(1,2);
end
