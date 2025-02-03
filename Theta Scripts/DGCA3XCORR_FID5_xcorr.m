%% Xcorr CA1 Subregion to G10 spikes FID 5
clear
clc
J6=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\times_J6.mat");
J6=J6.cluster_class(:,2);
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");
well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

DG_Elec=well_spike_dyn.channel_name(well_spike_dyn.regi==2 & well_spike_dyn.fi==6);

fs=25000;
t=1/fs:1/fs:300;
%% xcorr
xcorr_table=table();

J6_spike_train=zeros(1,length(t));
J6_spike_train(ismembertol(t,J6/1000))=1;

for nElec=1:length(DG_Elec)
    myElec=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33152 210715 21div 210806_1.h5\"+DG_Elec(nElec)+"_spikes.mat");
    my_spike_train=zeros(1,length(t));
    my_spike_train(ismembertol(t,myElec.index/1000))=1;
    [r,l]=xcorr(J6_spike_train,my_spike_train,fs*0.01,"normalized");
    figure
    plot(l/fs*1000,r)
    title(DG_Elec(nElec))

    [~,P]=corrcoef(J6_spike_train,my_spike_train);

    xcorr_table.elec(nElec)=DG_Elec(nElec);
    xcorr_table.r{nElec}=r;
    xcorr_table.l{nElec}=l;
    xcorr_table.pval(nElec)=P(1,2);
end
