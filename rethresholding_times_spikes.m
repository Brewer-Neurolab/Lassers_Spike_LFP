%% Spike rethresholding tunnels
clear
clc
close all
parent_dir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
myDir=dir(parent_dir);
folders=string({myDir.name});
folders=folders([myDir.isdir]);
folders=folders(3:end);

recording_time=300; %seconds
spike_direction="negative";

for fi=6%1:6

    subdir=dir(fullfile(parent_dir,folders(fi)));
    subdirFiles=string({subdir.name});
    timesFiles=subdirFiles(contains(subdirFiles,"times"));
    myElecs=split(timesFiles',[".","_"]);
    myElecs=myElecs(:,2);

    for nElec=22%1:length(myElecs)
        %load raw data
        rawData=load(fullfile(parent_dir,folders(fi),myElecs(nElec)+".mat"));
        fs=rawData.sr;
        t=0:1/fs:recording_time-(1/fs);
        rawData=rawData.data;
        %design wave clus filter
        [b,a]=ellip(2,0.1,40,[300,3000]/(fs/2));
        filtData=filtfilt(b,a,rawData);
        absData=abs(filtData);

        %load spike times
        mySpikes=load(fullfile(parent_dir,folders(fi),"times_"+myElecs(nElec)+".mat"));
        spike_times=mySpikes.cluster_class(:,2)/1000; %ms to s
        spike_idx=ismembertol(t,spike_times);
        spike_amplitudes=absData(spike_idx);

        %create threshold vector
        sigma=median(absData/0.6745);
        threshVec=(1:10)*sigma;

        %calculate the current spike amplitude distribution (from sigma=5)
        figure
        histogram(spike_amplitudes,logspace(log10(min(spike_amplitudes)),log10(max(spike_amplitudes)),40))
        set(gca,"XScale","log")
        title("FID "+fi+" "+myElecs(nElec))
        xline(threshVec)

        % first attempt ensure spikes are actually above threshold
        rethreshed_spikes=t(spike_idx & absData>=100*sigma)*1000;
        mySpikes.cluster_class=[ones(length(rethreshed_spikes),1),rethreshed_spikes'];

        %save in new dir
        if ~isfolder(fullfile(parent_dir,folders(fi),"rethreshed"))
            mkdir(fullfile(parent_dir,folders(fi),"rethreshed"))
        end
        
        save(fullfile(parent_dir,folders(fi),"rethreshed","times_"+myElecs(nElec)+".mat"),"-struct","mySpikes")
    end

end
