%% Read the data
clear all
close all
clc

eeglab

ControlGroup = [];

%ControlDir = 'e:\Dropbox (PUK)\Desktop\Microstate Demo\';
ControlDir = '/Users/Thomas/Desktop/Microstate Demo/';

Dir2Use = dir(fullfile(ControlDir,'*.vhdr'));

fnames = {Dir2Use.name};

for f = 1:numel(fnames)
    EEG = pop_fileio(fullfile(ControlDir,fnames{f}));
    EEG = eeg_RejectBABadIntervals( EEG);
    setname = strrep(fnames{f},'.vhdr','');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',fnames{f},'gui','off');
    EEG=pop_chanedit(EEG, 'lookup','/Users/Thomas/Documents/MATLAB/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
%    EEG=pop_chanedit(EEG, 'lookup','E:\Dropbox (PUK)\Docs\MATLAB\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
    EEG = pop_reref( EEG, []);
    EEG = pop_eegfiltnew(EEG, 2, 20, 424, 0, [], 0);
    ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET);
    ControlGroup = [ControlGroup CURRENTSET];
end

eeglab redraw
%% Cluster the stuff
ClustPars = struct('MinClasses',3,'MaxClasses',6,'GFPPeaks',true,'IgnorePolarity',true,'MaxMaps',500,'Restarts',20', 'UseAAHC',true);

for i = 1:numel(ControlGroup) 
    tmpEEG = eeg_retrieve(ALLEEG,ControlGroup(i));
    fprintf(1,'Clustering dataset %s (%i/%i)\n',tmpEEG.setname,i,numel(ControlGroup));
    tmpEEG = pop_FindMSTemplates(tmpEEG, ClustPars);
    ALLEEG = eeg_store(ALLEEG, tmpEEG, ControlGroup(i));
end

eeglab redraw

%% Now we combine across subjects and edit the mean
EEG = pop_CombMSTemplates(ALLEEG, ControlGroup, 0, 0, 'GrandMean Controls');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 
[ALLEEG,EEG] = pop_ShowIndMSMaps(EEG, 3, 1, ALLEEG);
[ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
GrandMeanControlIndex = CURRENTSET;

eeglab redraw
%% And we sort things out over subjects
ALLEEG = pop_SortMSTemplates(ALLEEG, ControlGroup, 0, GrandMeanControlIndex);
ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Save things
for f = 1:numel(ALLEEG)
    EEG = eeg_retrieve(ALLEEG,f);
    fname = [EEG.setname,'.vhdr'];
    EEG = pop_saveset( EEG, 'filename',fname,'filepath','E:\\Dropbox (PUK)\\Normative Database\\EEGLab stuff\\');
end

%% Visualize some stuff
FitPars = struct('nClasses',4,'lambda',0,'b',20,'PeakFit',false, 'BControl',true);
EEG = eeg_retrieve(ALLEEG,1); % First EEG...
pop_ShowIndMSDyn([],EEG,0,FitPars);
pop_ShowIndMSMaps(EEG,FitPars.nClasses);

%% And do the stats
FitPars = struct('nClasses',4,'lambda',0.3,'b',30,'PeakFit',true, 'BControl',false);

GrandMeanControlIndex = 474;
% Using the individual templates
pop_QuantMSTemplates(ALLEEG, ControlGroup, 0, FitPars, []                   , '/Users/Thomas/Desktop/test.csv');
%pop_QuantMSTemplates(ALLEEG, ControlGroup, 0, FitPars, []                   , 'E:\Dropbox (PUK)\Desktop\test.csv');
% Using the mean template from the controls
pop_QuantMSTemplates(ALLEEG, ControlGroup, 1, FitPars, GrandMeanControlIndex, 'E:\Dropbox (PUK)\Desktop\Test2.csv');

eeglab redraw
