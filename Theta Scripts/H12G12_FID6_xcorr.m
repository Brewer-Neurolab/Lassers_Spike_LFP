%% FID 6 H12 to G12 Crosscorrelation 10 ms

%load H12
H12=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\H12.mat");

%load G12
G12=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\G12.mat");

fs=25000;

% design filters
d=designfilt("bandpassiir",FilterOrder=6,SampleRate=fs,HalfPowerFrequency1=300,HalfPowerFrequency2=3000);
H12=filtfilt(d,H12.data);
G12=filtfilt(d,G12.data);

% xcorr all signal
[rAll,lAll]=xcorr(G12,H12,0.01*fs,"normalized");

%%
figure
plot(lAll/fs*1000,rAll)
xlabel("Lag (ms)")
ylabel("r")
xlim([-6,6])

max(rAll)
lAll(rAll==max(rAll))/fs

%% Xcorr for spike times
t=1/fs:1/fs:300;

%load H12
H12=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\4x 33168 210715 21div 210806_1.h5\H12_spikes.mat");
H12_spikes=H12.index;
H12_spike_train=zeros(1,length(t));
H12_idx=ismembertol(t,H12_spikes/1000);
H12_spike_train(H12_idx)=1;

%load G12
G12=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\times_G12.mat");
G12_spikes=G12.cluster_class(:,2);
G12_spike_train=zeros(1,length(t));
G12_idx=ismembertol(t,G12_spikes/1000);
G12_spike_train(G12_idx)=1;

[rSpikes,lSpikes]=xcorr(G12_spike_train,H12_spike_train,0.01*fs,"normalized");
[~,P]=corrcoef(G12_spike_train,H12_spike_train);

figure
plot(lSpikes/fs*1000,rSpikes,"LineWidth",1)
xlabel("Lag (ms)")
ylabel("r")
xlim([-6,6])
set(gca,"FontSize",16)