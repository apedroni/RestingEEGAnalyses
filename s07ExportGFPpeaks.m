%% Export GFPpeaks for Cartool
clear, clc;

warning('off','all')
% this is the path with the main analyses
workingDirectory = pwd;

% this is the path with the results:
rootpath =  '/Volumes/methlab/Neurometric/Test_Retest_Data/';
% path eye-tracker files
% add EEGLAB path
addpath('~/Dropbox/EEG_analysis/GeneralMatlab/eeglab14_1_1b/')
eeglab
close
% add functions paths
addpath('~/Dropbox/AA_Neurometric/ANALYSES/functions/Rest/');
addpath('~/Dropbox/AA_Neurometric/ANALYSES/functions/')
addpath('~/Dropbox/AA_Neurometric/ANALYSES/')

%% get the complete Datasets (with two timepoints) 
load('~/Dropbox/AA_Neurometric/ANALYSES/RetestResults/RestingEEGwithMicrostates.mat')

i = 1;
GFPpeaksT1 = []; GFPpeaksT2 = [];
for k=1:length(AllData)
    if ~isempty(AllData(k).T1) &&  ~isempty(AllData(k).T2)...
            && ~strcmp(AllData(k).T1.ID(1),'b') && ~strcmp(AllData(k).T2.ID(1),'b') %% and if its a bad 
         CompleteData(i).T1 = AllData(k).T1;
         CompleteData(i).T2 = AllData(k).T2;
        
        GFPpeaksT1 = [GFPpeaksT1   AllData(k).T1.microstate.data];
        GFPpeaksT2 = [GFPpeaksT2   AllData(k).T2.microstate.data];
        i = i + 1;
    end
end


dlmwrite(['~/Dropbox/AA_Neurometric/ANALYSES/RestingEEGAnalyses/MicrostateResults/','GFPpeaksT1.ep'], GFPpeaksT1','delimiter',' ')
dlmwrite(['~/Dropbox/AA_Neurometric/ANALYSES/RestingEEGAnalyses/MicrostateResults/','GFPpeaksT2.ep'], GFPpeaksT2','delimiter',' ')




