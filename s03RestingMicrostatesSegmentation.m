%% Microstates Analysis 
% use GEEG with X GFP peak maps as input for segmentation
clear, clc;

RestingCreateSettings;
%% Backfitting
workingDirectory = '~/Dropbox/AA_Neurometric/ANALYSES/RestingEEGAnalyses/';

% this is the path with the data:
eegpath =  '/Volumes/methlab/Neurometric/Test_Retest_Data_Neurometric_Sahar_New_results/';


% add EEGLAB path
addpath('~/Dropbox/EEG_analysis/GeneralMatlab/eeglab14_1_1b/')
eeglab
close

% add functions paths
addpath('~/Dropbox/AA_Neurometric/ANALYSES/functions/Rest/');
addpath('~/Dropbox/AA_Neurometric/ANALYSES/functions/')
addpath('~/Dropbox/AA_Neurometric/ANALYSES/')

load('~/Dropbox/AA_Neurometric/ANALYSES/RetestResults/chanlocs105.mat')
load('~/Dropbox/AA_Neurometric/ANALYSES/RetestResults/RestingEEG.mat')

i = 1;
for k=1:length(AllData)
    if ~isempty(AllData(k).T1) &&  ~isempty(AllData(k).T2) ...
            && ~strcmp(AllData(k).T1.ID(1),'b') && ~strcmp(AllData(k).T2.ID(1),'b') %% and if its a bad 

        CompleteData(i).T1 = AllData(k).T1;
        CompleteData(i).T2 = AllData(k).T2;      
        i = i + 1;
           
    end
end
%% A) SEGMENTATION on concatenated GFP peak maps OUR STYLE

%% T1 
EEGtmp = eeg_emptyset();
EEGtmp.setname = 'GFPpeakmaps'
EEGtmp.chanlocs = chanlocs;
EEGtmp.nbchan = length(chanlocs)
EEGtmp.trials = 1
EEGtmp.srate = 500;

for i = 1:length(CompleteData) 
EEGtmp.data = [EEGtmp.data , CompleteData(i).T1.microstate.data];
end
EEGtmp.microstate.data = EEGtmp.data;
EEG.pnts = size(EEGtmp.data,2)
EEGtmp.times = (1:size(EEG.data,2))*1000/EEG.srate

%% 2. Segment Data into Microstates
EEG = pop_micro_segment ( EEGtmp, 'algorithm','modkmeans', 'Nmicrostates', 4,'Nrepetitions',50)

h = figure
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(EEG.microstate.prototypes(:,z)',chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsT1.png');

T1 = EEG.microstate;

%% T2 
EEG = eeg_emptyset();
EEG.setname = 'GFPpeakmaps'
EEG.chanlocs = chanlocs;
EEG.nbchan = length(chanlocs)
EEG.trials = 1
EEG.srate = 500;

for i = 1:length(CompleteData) 
EEG.data = [EEG.data , CompleteData(i).T2.microstate.data];
end

EEG.microstate.data = EEG.data;
EEG.pnts = size(EEG.data,2);
EEG.times = (1:size(EEG.data,2))*1000/EEG.srate;

%% 2. Segment Data into Microstates
EEG = pop_micro_segment ( EEG, 'algorithm','modkmeans', 'Nmicrostates', 4,'Nrepetitions',50)

h = figure
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(EEG.microstate.prototypes(:,z)',chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsT2.png');

T2 = EEG.microstate



%% T1 and T2 combined

EEG = eeg_emptyset();
EEG.setname = 'GFPpeakmaps'
EEG.chanlocs = chanlocs;
EEG.nbchan = length(chanlocs)
EEG.trials = 1
EEG.srate = 500;

for i = 1:length(CompleteData) 
EEG.data = [EEG.data ,  CompleteData(i).T1.microstate.data, CompleteData(i).T2.microstate.data ];
end
EEG.microstate.data =  EEG.data;
EEG.pnts = size(EEG.data,2)
EEG.times = (1:size(EEG.data,2))*1000/EEG.srate

%% 2. Segment Data into Microstates
EEG = pop_micro_segment ( EEG, 'algorithm','modkmeans', 'Nmicrostates', 4,'Nrepetitions',50)

h = figure
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(EEG.microstate.prototypes(:,z)',chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsT1T2.png');

T1T2 = EEG.microstate

Prototypes.chanlocs = chanlocs;
Prototypes.OurVersion.T1 = T1;
Prototypes.OurVersion.T2 = T2;
Prototypes.OurVersion.T1T2 = T1T2;


%% B) Sorting the microstates with thomas' algorithm... 
%load(['~/Dropbox/AA_Neurometric/ANALYSES/GroupLevelResults/Resting/' 'SingleSubjectSegmentation.mat'])

SingleProtos = []
for i = 1:length(CompleteData) 
SingleProtos(i,:,:) = CompleteData(i).T1.microstate.prototypes';
end

tmp2 = squeeze(SingleProtos(2,1,:))

[MeanMap,SortedMaps,OldMapFit] = PermutedMeanMaps(double(SingleProtos),0);

h = figure;
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(MeanMap(z,:),chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsThomasVersionT1.png');

% plot all Microstate Segmentations
% figure;
% o = 1;
% for k = 1:size(SortedMaps,1)
%     for z = 1:4
%         map = squeeze(SortedMaps(k,:,:));
%         subplot(size(SortedMaps,1),4,o)
%         topoplot(map(z,:)',chanlocs ,'style','fill','conv','on','electrodes','off')
%         o = o + 1 ;
%     end
% end
% for z = 1:4
%     subplot(size(SortedMaps,1)+1,4,z+o-1)
%     title(['Map ' num2str(z)])
%     topoplot(MeanMap(z,:)',chanlocs,'style','fill','conv','on','electrodes','off' )
% end

for i = 1:length(CompleteData)
Subjects{i,:} = CompleteData(i).T1.ID
end


Prototypes.Thomas.T1.prototypes =  MeanMap;
Prototypes.Thomas.T1.sortedprototypes =  SortedMaps;
Prototypes.Thomas.T1.Subjects = Subjects

%save('./MicrostateResults/ProtosSubThomasT1.mat','MeanMap','SortedMaps','-v7.3');
   


SingleProtos = []
for i = 1:length(CompleteData) 
SingleProtos(i,:,:) = CompleteData(i).T2.microstate.prototypes';
end

for i = 1:length(CompleteData)
Subjects{i,:} = CompleteData(i).T2.ID
end

[MeanMap,SortedMaps,OldMapFit] = PermutedMeanMaps(double(SingleProtos),0);

h = figure;
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(MeanMap(z,:),chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsThomasVersionT2.png');

% plot all Microstate Segmentations
% figure;
% o = 1;
% for k = 1:size(SortedMaps,1)
%     for z = 1:4
%         map = squeeze(SortedMaps(k,:,:));
%         subplot(size(SortedMaps,1),4,o)
%         topoplot(map(z,:)',chanlocs ,'style','fill','conv','on','electrodes','off')
%         o = o + 1 ;
%     end
% end
% for z = 1:4
%     subplot(size(SortedMaps,1)+1,4,z+o-1)
%     title(['Map ' num2str(z)])
%     topoplot(MeanMap(z,:)',chanlocs,'style','fill','conv','on','electrodes','off' )
% end

for i = 1:length(CompleteData)
Subjects{i,:} = CompleteData(i).T2.ID
end

Prototypes.Thomas.T2.prototypes =  MeanMap;
Prototypes.Thomas.T2.sortedprototypes =  SortedMaps;
Prototypes.Thomas.T2.Subjects = Subjects

%save('./MicrostateResults/ProtosSubThomasT2.mat','MeanMap','SortedMaps','-v7.3');


%% T1T2


for i = 1:length(CompleteData) 
SingleProtos1(i,:,:) = CompleteData(i).T1.microstate.prototypes';
SingleProtos2(i,:,:) = CompleteData(i).T2.microstate.prototypes';
end

SingleProtos = [SingleProtos1 ; SingleProtos2]

[MeanMap,SortedMaps,OldMapFit] = PermutedMeanMaps(double(SingleProtos),0);

h = figure;
for z = 1:4
    subplot(1,4,z)
    title(['Map ' num2str(z)])
    topoplot(MeanMap(z,:),chanlocs )
end
saveas(h,'./MicrostateResults/ProtoMapsThomasVersionT2.png');

% plot all Microstate Segmentations
% figure;
% o = 1;
% for k = 1:size(SortedMaps,1)
%     for z = 1:4
%         map = squeeze(SortedMaps(k,:,:));
%         subplot(size(SortedMaps,1),4,o)
%         topoplot(map(z,:)',chanlocs ,'style','fill','conv','on','electrodes','off')
%         o = o + 1 ;
%     end
% end
% for z = 1:4
%     subplot(size(SortedMaps,1)+1,4,z+o-1)
%     title(['Map ' num2str(z)])
%     topoplot(MeanMap(z,:)',chanlocs,'style','fill','conv','on','electrodes','off' )
% end

Prototypes.Thomas.T1T2.prototypes =  MeanMap;
Prototypes.Thomas.T1T2.sortedprototypes =  SortedMaps;

%save('./MicrostateResults/ProtosSubThomasT1T2.mat','MeanMap','SortedMaps','-v7.3');
   

save('~/Dropbox/AA_Neurometric/ANALYSES/RestingEEGAnalyses/MicrostateResults/Prototypes.mat','Prototypes','-v7.3') 
