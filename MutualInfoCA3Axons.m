%% SETUP

clear
clc
close all

parentDir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
dirList=dir(parentDir);
subDirs=string({dirList.name});
subDirs=subDirs([dirList.isdir]==1);
subDirs=subDirs(3:end);

axonSpikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axonSpikes=axonSpikes.allregion_unit_matched_stim;

fs=25000;
re_fs=1000;
t_rec=300;
%% Cycle through FF axons and make heat map (NO FB)

for FI=1:length(subDirs)
    MIMat{FI}=NaN([4,4]);
    CA3Rows=find(axonSpikes{FI}.Subregion=="CA3-CA1");
    %check if ff
    for n1=1:length(CA3Rows)
        for n2=1:length(CA3Rows)
            if ~isempty(axonSpikes{FI}.up_ff(CA3Rows(n1))) & ~isempty(axonSpikes{FI}.up_ff(CA3Rows(n2))) & n1~=n2
                axon1=strsplit(axonSpikes{FI}.("Electrode Pairs")(CA3Rows(n1)),{'-'});
                axon1Dir=fullfile(parentDir,subDirs(FI),axon1(1)+".mat");
                axon2=strsplit(axonSpikes{FI}.("Electrode Pairs")(CA3Rows(n2)),{'-'});
                axon2Dir=fullfile(parentDir,subDirs(FI),axon2(1)+".mat");
                MIMat{FI}(n1,n2)=twoAxonMutInfo(axon1Dir,axon2Dir,4,10,fs,re_fs,t_rec,100);
            end
        end
    end
    figure
    xvals=split(axonSpikes{FI}.("Electrode Pairs")(CA3Rows(1:4)),{'-'});
    xvals=xvals(:,1);
    yvals=xvals;
    heatmap(xvals,yvals,round(MIMat{FI},2))
    title("Mutual Info (Bits) FID "+FI)
    set(gca,"FontSize",16)

    disp("FID "+FI)
end
%% Histogram of MI MAT

HistVec=reshape(cell2mat(MIMat),[],1);
figure
histogram(HistVec)
xlabel("Mutual Information (Bits)")
ylabel("Axon Pair Counts")
set(gca,"FontSize",16)