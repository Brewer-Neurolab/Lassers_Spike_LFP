%% SETUP

clear
clc
close all

saveDir="C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Images\Theta\MutualInfoHeatmaps";

parentDir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
dirList=dir(parentDir);
subDirs=string({dirList.name});
subDirs=subDirs([dirList.isdir]==1);
subDirs=subDirs(3:end);

axonSpikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axonSpikes=axonSpikes.allregion_unit_matched_stim;

interRegions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];
subregions=["EC","DG","CA3","CA1"];

% for fi=1:length(axonSpikes)
%     for nelec=1:height(axonSpikes{fi})
%         if ~isempty(axonSpikes{fi}.up_ff{nelec}) | (~isempty(axonSpikes{fi}.up_fb{nelec}) & axonSpikes{fi}.Subregion(nelec)=="CA1-EC")
%             if isempty(axonSpikes{fi}.up_ff{nelec}) & ~isempty(axonSpikes{fi}.up_fb{nelec}) & axonSpikes{fi}.Subregion(nelec)=="CA1-EC"
%                 axonSpikes{fi}.Subregion(nelec)="EC-CA1";
%             end
%         end
%     end
% end

%3.5 SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

all_reg=[interRegions];

fs=25000;
re_fs=1000;
t_rec=300;
%% Cycle through FF axons and make heat map (NO FB)

for nReg=1:length(all_reg)
    clear MIMat
    for FI=1:length(subDirs)
        MIMat{FI}=NaN([4,4]);
        Rows=find(axonSpikes{FI}.Subregion==all_reg(nReg));
        if axonSpikes{FI}.Subregion=="EC-CA3"
            MIMat{FI}=NaN([5,5]);
        end
        %check if ff
        for n1=1:length(Rows)
            for n2=1:length(Rows)
                if n2<=n1
                    MIMat{FI}(n1,n2)=NaN;
                    continue
                end

                % if (~isempty(axonSpikes{FI}.up_ff{Rows(n1)}) & ~isempty(axonSpikes{FI}.up_ff{Rows(n2)})) | (~isempty(axonSpikes{FI}.up_fb{Rows(n1)}) & ~isempty(axonSpikes{FI}.up_fb{Rows(n2)}))
                    axon1=strsplit(axonSpikes{FI}.("Electrode Pairs")(Rows(n1)),{'-'});
                    axon1Dir=fullfile(parentDir,subDirs(FI),axon1(1)+".mat");
                    axon2=strsplit(axonSpikes{FI}.("Electrode Pairs")(Rows(n2)),{'-'});
                    axon2Dir=fullfile(parentDir,subDirs(FI),axon2(1)+".mat");
                    MIMat{FI}(n1,n2)=twoAxonMutInfo(axon1Dir,axon2Dir,4,10,fs,re_fs,t_rec,100);
                % else
                    % MIMat{FI}(n1,n2)=NaN;
                % end
            end
        end
        figure('units','normalized','outerposition',[0 0 1 1])
        xvals=split(axonSpikes{FI}.("Electrode Pairs")(Rows),{'-'});
        xvals=xvals(:,1);
        yvals=xvals;
        k=heatmap(xvals,yvals,round(MIMat{FI},2));
        k.Position=[.25 .15 .5 .8];
        % title("Mutual Info (Bits) FID "+FI)
        set(gca,"FontSize",50)

        disp("FID "+FI)
        
        saveas(k,fullfile(saveDir,all_reg(nReg)+" FID "+FI),"png")
        close all
    end
    
end
disp("FF Calc Done!")
%% EC-CA1 FB

% SETUP

clear
clc
close all

saveDir="C:\Users\lasss\Documents\Research\Brewer Lab work\Code\Lassers_Spike_LFP\Images\Theta\MutualInfoHeatmaps";

parentDir="D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A";
dirList=dir(parentDir);
subDirs=string({dirList.name});
subDirs=subDirs([dirList.isdir]==1);
subDirs=subDirs(3:end);

axonSpikes=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\18-Apr-2023_A\allregion_unit_matched_cleaned.mat");
axonSpikes=axonSpikes.allregion_unit_matched_stim;

interRegions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC","EC-CA3"];
subregions=["EC","DG","CA3","CA1"];

for fi=1:length(axonSpikes)
    for nelec=1:height(axonSpikes{fi})
        if axonSpikes{fi}.Subregion(nelec)=="CA1-EC"
            axonSpikes{fi}.Subregion(nelec)="EC-CA1";
        end
    end
end

%3.5 SD
% well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes\well_spike_dynamics_table_hfs_3-5.mat");

%5SD
well_spike_dyn=load("D:\Brewer lab data\Slow Oscillation 4 Chamber 5 Tunnel Arrays\4x 210715 210806\1\Well Spikes 5SD Min\well_spike_dynamics_table_hfs.mat");

well_spike_dyn=well_spike_dyn.well_spike_dynamics_table;

all_reg=[interRegions];

fs=25000;
re_fs=1000;
t_rec=300;

%% Cycle through FB in EC-CA1 axons and make heat map 

myReg="EC-CA1";

for nReg=1:length(myReg)
    clear MIMat
    for FI=1:length(subDirs)
        MIMat{FI}=NaN([4,4]);
        Rows=find(axonSpikes{FI}.Subregion==myReg(nReg));
        if axonSpikes{FI}.Subregion=="EC-CA3"
            MIMat{FI}=NaN([5,5]);
        end
        %check if ff
        for n1=1:length(Rows)
            for n2=1:length(Rows)
                if n2<=n1
                    MIMat{FI}(n1,n2)=NaN;
                    continue
                end
                % if ~isempty(axonSpikes{FI}.up_fb{Rows(n1)}) & ~isempty(axonSpikes{FI}.up_fb{Rows(n2)}) & n1~=n2
                    axon1=strsplit(axonSpikes{FI}.("Electrode Pairs")(Rows(n1)),{'-'});
                    axon1Dir=fullfile(parentDir,subDirs(FI),axon1(1)+".mat");
                    axon2=strsplit(axonSpikes{FI}.("Electrode Pairs")(Rows(n2)),{'-'});
                    axon2Dir=fullfile(parentDir,subDirs(FI),axon2(1)+".mat");
                    MIMat{FI}(n1,n2)=twoAxonMutInfo(axon1Dir,axon2Dir,4,10,fs,re_fs,t_rec,100);
                % else
                %     MIMat{FI}(n1,n2)=NaN;
                % end
            end
        end
        f=figure;
        xvals=split(axonSpikes{FI}.("Electrode Pairs")(Rows),{'-'});
        xvals=xvals(:,1);
        yvals=xvals;
        heatmap(xvals,yvals,round(MIMat{FI},2))
        title("Mutual Info (Bits) FID "+FI)
        set(gca,"FontSize",16)

        disp("FID "+FI)
        saveas(f,fullfile(saveDir,myReg(nReg)+" FID "+FI),"png")
        close all
    end

end
disp("FB Calc Done!")