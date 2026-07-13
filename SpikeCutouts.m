%% Copy and paste a spike at theta frequency 

clear
clc
close all

% data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33152 210715 21div 210806_1_mat_files\G12.mat");
data=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\4x 33168 210715 21div 210806_1_mat_files\G12.mat");
fs=data.sr;
data=data.data;

[A,B,C,D]=butter(8,3000/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
spike_data=filtfilt(sos,g,data);

[A,B,C,D]=butter(8,300/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
spike_data=filtfilt(sos,g,spike_data);

t=0:1/fs:300-(1/fs);

% [~,spike_bound_L]=min(abs(t-11.514));
% [~,spike_bound_R]=min(abs(t-11.5166));
[~,spike_bound_L]=min(abs(t-15.2186));
[~,spike_bound_R]=min(abs(t-15.2208));

spike_cutout=spike_data(spike_bound_L:spike_bound_R);
t_cutout=1/fs:1/fs:length(spike_cutout)/fs;

spike_cutout_padded=[zeros(1,1*fs),spike_cutout,zeros(1,1*fs)];
t_padded=1/fs:1/fs:length(spike_cutout_padded)/fs;


figure
plot(spike_cutout_padded)

% spectrogram(spike_cutout_padded,[],[],100,fs,"yaxis");
% [cfs,frq]=cwt(spike_cutout_padded,"morse",fs);
% fb=cwtfilterbank(SignalLength=numel(spike_cutout_padded),SamplingFrequency=fs,FrequencyLimits=[1 300]);
% cwt(spike_cutout_padded,"ExtendSignal",1,"Frequency",fs,"FrequencyLimits",[0.5,300])
[cfs,frq]=cwt(spike_cutout_padded,"amor",fs,"FrequencyLimits",[0.5,3000]);

figure
surface(t_padded,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,4,10,30,100,300,1000,3000])
%% Raw signal spectrum
% sample_raw=data(270000:320000);
sample_raw=data(250000:318500);
[cfs,frq]=cwt(sample_raw,"amor",fs,"FrequencyLimits",[3,300]);
t_rep=1/fs:1/fs:length(sample_raw)/fs;
figure
plot(t_rep,sample_raw,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")
xlim([min(t_rep),max(t_rep)])
ax=gca;
set(ax,"FontSize",50)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
xlim([0.5,2])
axis square


figure
surface(t_rep,frq,log10(abs(cfs)./frq))
% surface(t_rep,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([3,10,30,100,300])
c=colorbar;
c.Label.String="Log10 Magnitude/Hz";
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])
ylim([3,300])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"FontSize",50)
% colormap("jet")
colormap(flipud(slanCM("plasma")))
axis square

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_sample=filtfilt(sos,g,sample_raw);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_sample=filtfilt(sos,g,filtered_sample);

figure
plot(t_rep,filtered_sample,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")

hilbertAmp=abs(hilbert(filtered_sample));
hold on
plot(t_rep,hilbertAmp,'Color','g','LineWidth',2)
hold off
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])
ax=gca;
set(ax,"FontSize",50)
% xticks(-1000:500:1000)
% xticklabels(-1000:500:1000)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
axis square

% thresh=std(hilbertAmp);
% aboveThresh=find(hilbertAmp>thresh);
thetaPowerCFS=abs(cfs(frq>=4 & frq<=10,:))./frq(frq>=4 & frq<=10,:);
% thetaPowerThreshed=thetaPowerCFS(:,aboveThresh(1):aboveThresh(end));

thetaPowerRaw=sum(thetaPowerCFS,"all");
%% Spike repeats 3

spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
t_rep=1/fs:1/fs:length(spike_repeat)/fs;
figure
plot(spike_repeat)

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[0.5,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
surface(t_rep,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,4,10,30,100,300])

%% Spike repeats 10 at 5Hz
% spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
% t_rep=1/fs:1/fs:length(spike_repeat)/fs;
[spike_repeat,t_rep]=spike_repeater(spike_cutout,0.2,fs,10,0,1);
figure
plot(t_rep,spike_repeat)
xlabel("Time (s)")
ylabel("uV")

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[0.5,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,10*log10(abs(cfs)))
% surface(t_rep,frq,10*log10(abs(cfs))./repmat(frq,[1,size(cfs,2)]))
surface(t_rep,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,4,10,30,100,300])
xlabel("Time (s)")
ylabel("Frequency")
% colormap("jet")
colormap(flipud(slanCM("plasma")))
colorbar

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,spike_repeat);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,filtered_spikes);

figure
plot(t_rep,filtered_spikes)
xlabel("Time (s)")
ylabel("uV")

hilbertAmp=abs(hilbert(filtered_spikes));
hold on
plot(t_rep,hilbertAmp)
hold off
%% Spike repeats 3 at 5Hz
% spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
% t_rep=1/fs:1/fs:length(spike_repeat)/fs;
[spike_repeat,t_rep]=spike_repeater(spike_cutout,0.2,fs,3,0,1);
figure
plot(t_rep,spike_repeat,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")
ax=gca;
set(ax,"FontSize",50)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
xlim([0.5,2])
axis square

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[3,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,10*log10(abs(cfs)))
% surface(t_rep,frq,10*log10(abs(cfs))./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,abs(cfs)./frq)
surface(t_rep,frq,log10(abs(cfs)./frq))
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,3,10,30,100,300])
xlabel("Time (s)")
ylabel("Frequency")
% colormap("jet")
colormap(flipud(slanCM("plasma")))
c=colorbar;
clim([-4,2])
set(gca,"FontSize",50)
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])
ylim([3,300])
c.Label.String="Log Magnitude/Hz";
axis square

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,spike_repeat);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,filtered_spikes);

figure
plot(t_rep,filtered_spikes,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])

hilbertAmp=abs(hilbert(filtered_spikes));
hold on
plot(t_rep,hilbertAmp,'Color','g','LineWidth',2)
hold off
ax=gca;
set(ax,"FontSize",50)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
axis square

thresh=std(hilbertAmp);
% aboveThresh=find(hilbertAmp>thresh);
thetaPowerCFS=abs(cfs(frq>=4 & frq<=10,:))./frq(frq>=4 & frq<=10,:);
% thetaPowerThreshed=thetaPowerCFS(:,aboveThresh(1):aboveThresh(end));

thetaPower5Hz=sum(thetaPowerCFS,"all");
%% Spike repeats 22 at 200 Hz
% spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
% t_rep=1/fs:1/fs:length(spike_repeat)/fs;
[spike_repeat,t_rep]=spike_repeater(spike_cutout,0.005,fs,22,0,1);
figure
plot(t_rep,spike_repeat,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")
xlim([min(t_rep),max(t_rep)])
ax=gca;
set(ax,"FontSize",18)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[0.5,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,10*log10(abs(cfs)))
% surface(t_rep,frq,10*log10(abs(cfs))./repmat(frq,[1,size(cfs,2)]))
surface(t_rep,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,4,10,30,100,300])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
% colormap("jet")
colormap(flipud(slanCM("plasma")))
c=colorbar;
c.Label.String="Magnitude/Hz";
xlim([min(t_rep),max(t_rep)])
ylim([min(frq),max(frq)])
set(gca,"FontSize",18)

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,spike_repeat);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,filtered_spikes);

figure
plot(t_rep,filtered_spikes,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")

hilbertAmp=abs(hilbert(filtered_spikes));
hold on
plot(t_rep,hilbertAmp,'Color','g','LineWidth',2)
hold off
xlim([min(t_rep),max(t_rep)])
ax=gca;
set(ax,"FontSize",18)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";

thresh=std(hilbertAmp);
% aboveThresh=find(hilbertAmp>thresh);
thetaPowerCFS=abs(cfs(frq>=4 & frq<=10,:))./frq(frq>=4 & frq<=10,:);
% thetaPowerThreshed=thetaPowerCFS(:,aboveThresh(1):aboveThresh(end));

thetaPower200Hz=sum(thetaPowerCFS,"all");
% yline(thresh)
%% Spike repeats 5 at 200 Hz
% spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
% t_rep=1/fs:1/fs:length(spike_repeat)/fs;
[spike_repeat,t_rep]=spike_repeater(spike_cutout,0.005,fs,5,0.2,3);
figure
plot(t_rep,spike_repeat,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")
xlim([min(t_rep),max(t_rep)])
ax=gca;
set(ax,"FontSize",50)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
xlim([0.5,2])
axis square

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[3,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,10*log10(abs(cfs)))
% surface(t_rep,frq,10*log10(abs(cfs))./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,abs(cfs)./frq)
surface(t_rep,frq,log10(abs(cfs)./frq))
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,3,10,30,100,300])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
% colormap("jet")
colormap(flipud(slanCM("plasma")))
c=colorbar;
clim([-4,2])
c.Label.String="Log Magnitude/Hz";
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])
ylim([3,300])
set(gca,"FontSize",50)
% colormap("jet")
colormap(flipud(slanCM("plasma")))
axis square

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,spike_repeat);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,filtered_spikes);

figure
plot(t_rep,filtered_spikes,'Color','k','LineWidth',2)
xlabel("Time (s)")
ylabel("uV")

hilbertAmp=abs(hilbert(filtered_spikes));
hold on
plot(t_rep,hilbertAmp,'Color','g','LineWidth',2)
hold off
% xlim([min(t_rep),max(t_rep)])
xlim([0.5,2])
ax=gca;
set(ax,"FontSize",50)
ax.TickLength(1)=0.025;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
axis square

thresh=std(hilbertAmp);
% aboveThresh=find(hilbertAmp>thresh);
thetaPowerCFS=abs(cfs(frq>=4 & frq<=10,:))./frq(frq>=4 & frq<=10,:);
% thetaPowerThreshed=thetaPowerCFS(:,aboveThresh(1):aboveThresh(end));

thetaPower200Hz=sum(thetaPowerCFS,"all");
% yline(thresh)
%% Spike repeats 1 
% spike_repeat=[zeros(1,1*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,0.2*fs),spike_cutout,zeros(1,1*fs)];
% t_rep=1/fs:1/fs:length(spike_repeat)/fs;
[spike_repeat,t_rep]=spike_repeater(spike_cutout,0.005,fs,1,0,1);
figure
plot(t_rep,spike_repeat)
xlabel("Time (s)")
ylabel("uV")
xlim([min(t_rep),max(t_rep)])
set(gca,"FontSize",18)

[cfs,frq]=cwt(spike_repeat,"amor",fs,"FrequencyLimits",[0.5,300]);

figure
% surface(t_rep,frq,abs(cfs)./repmat(frq,[1,size(cfs,2)]))
% surface(t_rep,frq,10*log10(abs(cfs)))
% surface(t_rep,frq,10*log10(abs(cfs))./repmat(frq,[1,size(cfs,2)]))
surface(t_rep,frq,abs(cfs)./frq)
shading flat
ax=gca;
ax.YScale="log";
yticks([0.5,4,10,30,100,300])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
% colormap("jet")
colormap(flipud(slanCM("plasma")))
c=colorbar;
c.Label.String="Magnitude/Hz";
xlim([min(t_rep),max(t_rep)])
set(gca,"FontSize",18)

%power from hilbert
[A,B,C,D]=butter(8,10/(fs/2),'low');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,spike_repeat);

[A,B,C,D]=butter(8,4/(fs/2),'high');
[sos,g]=ss2sos(A,B,C,D);
filtered_spikes=filtfilt(sos,g,filtered_spikes);

figure
plot(t_rep,filtered_spikes)
xlabel("Time (s)")
ylabel("uV")

hilbertAmp=abs(hilbert(filtered_spikes));
hold on
plot(t_rep,hilbertAmp)
hold off
xlim([min(t_rep),max(t_rep)])
set(gca,"FontSize",18)

thresh=std(hilbertAmp);
% aboveThresh=find(hilbertAmp>thresh);
thetaPowerCFS=abs(cfs(frq>=4 & frq<=10,:))./frq(frq>=4 & frq<=10,:);
% thetaPowerThreshed=thetaPowerCFS(:,aboveThresh(1):aboveThresh(end));

thetaPower1=sum(thetaPowerCFS,"all");
% yline(thresh)
%% Bar chart of theta powers
% Powers are divided by frequency for normalization
figure
cats=categorical(["Raw","Single","5 Hz","200 Hz x3"]);
cats=reordercats(cats,["Raw","Single","5 Hz","200 Hz x3"]);
bar(cats,[thetaPowerRaw,thetaPower1,thetaPower5Hz,thetaPower200Hz])
ax=gca;
ax.YScale="log";
ylim([10^3,10^8])
ax.FontSize=14;
ylabel("Power uV^2")
ax=gca;
set(ax,"FontSize",24)
ax.TickLength(1)=0.05;
ax.LineWidth=2;
ax.XMinorTick="on";
ax.YMinorTick="on";
%% Functions
function [burst,t_rep]=spike_repeater(spike_cutout,space,fs,n,train_space,nTrain)

burst=zeros(1,1*fs);

for trains=1:nTrain
    for i=1:n
        if i==n
            burst=[burst,spike_cutout];
        else
            burst=[burst,spike_cutout,zeros(1,space*fs)];
        end
    end
    burst=[burst,zeros(1,train_space*fs)];
end
burst=[burst,zeros(1,1*fs)];

t_rep=1/fs:1/fs:length(burst)/fs;
end