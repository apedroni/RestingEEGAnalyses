%% Demo script for microstate analyses in EEGLAB
%
% %Author: Thomas Koenig, University of Bern, Switzerland, 2017
%  
%   Copyright (C) 2017 Thomas Koenig, University of Bern, Switzerland, 2016
%   thomas.koenig@puk.unibe.ch
%  
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%  
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%  
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
% ---------------------------------

% This is a sample script that you will have to adapt to meet your specific
% needs.

%% Define the basic parameters
% This is for vision analyzer data and may need adjustments
clear all
close all
clc

LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;

% For already saved and filtered EEG-lab data
% ReadVision = false;
% FilterTheData = false;

% for "fresh" vision analyzer data:
ReadVision = false;
FilterTheData = false;


% These are the paramters for the fitting based on GFP peaks only
FitPars = struct('nClasses',4,'lambda',1,'b',20,'PeakFit',false, 'BControl',true,'Rectify',false,'Normalize',false);

% Define the parameters for clustering
ClustPars = struct('MinClasses',4,'MaxClasses',4,'GFPPeaks',true,'IgnorePolarity',true,'MaxMaps',500,'Restarts',20', 'UseAAHC',false);
 
nGroups = str2double(inputdlg('Number of groups','Microstate analysis',1));

%% Read the data

eeglabpath = fileparts(which('eeglab.m'));                  % Getting eeglab path
DipFitPath = fullfile(eeglabpath,'plugins','dipfit2.3');  % Crafting DipFit Path

SavePath   = uigetdir([],'Path to store the results')

[SaveFileName,SavePathName] = uiputfile('*.mat','Save everything to one file')

if SavePath == 0
    return
end

eeglab

AllSubjects = [];
for Group = 1:nGroups

    GroupDir = uigetdir([],sprintf('Path to the data of group %i (Vision Analyzer data)',Group))
    
    if GroupDir == 0
        return
    end
    
    GroupIndex{Group} = [];
    
    if ReadVision == true
        DirGroup = dir(fullfile(GroupDir,'*.vhdr'));
    else
        DirGroup = dir(fullfile(GroupDir,'*.set'));
    end

    FileNamesGroup = {DirGroup.name};

    % Read the data from the group 
    for f = 1:numel(FileNamesGroup)
        if ReadVision == true
            tmpEEG = pop_fileio(fullfile(GroupDir,FileNamesGroup{f}));   % Basic file read
            tmpEEG = eeg_RejectBABadIntervals(tmpEEG);   % Get rid of bad intervals
            setname = strrep(FileNamesGroup{f},'.vhdr',''); % Set a useful name of the dataset
            [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'setname',FileNamesGroup{f},'gui','off'); % And make this a new set
            tmpEEG=pop_chanedit(tmpEEG, 'lookup',fullfile(DipFitPath,'standard_BESA','standard-10-5-cap385.elp')); % Add the channel positions
        else
            tmpEEG = pop_loadset('filename',FileNamesGroup{f},'filepath',GroupDir);
            [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off','setname',FileNamesGroup{f}) % And make this a new set
        end

        tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
        if FilterTheData == true
            tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % And bandpass-filter 2-20Hz
        end
        tmpEEG.group = sprintf('Group_%i',Group); % Set the group (will appear in the statistics output)
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
        GroupIndex{Group} = [GroupIndex{Group} CURRENTSET]; % And keep track of the groups
        AllSubjects = [AllSubjects CURRENTSET];
    end
    
end

eeglab redraw

save(fullfile(SavePathName,SaveFileName));
   
%% Cluster the stuff

load('~/Dropbox/AA_Neurometric/ANALYSES/RestingEEGAnalyses/ABIMAnalyses/Thomas/SegmentationThomas.mat')
ThomasVersion = EEG.msinfo.MSMaps(4).Maps';


% Loop across all subjects to identify the individual clusters
% for i = 1:numel(AllSubjects ) 
    tmpEEG = eeg_retrieve(ALLEEG,AllSubjects(1)); % the EEG we want to work with
    %fprintf(1,'Clustering dataset %s (%i/%i)\n',tmpEEG.setname,i,numel(AllSubjects )); % Some info for the impatient user
    tmpEEG = pop_FindMSTemplates(tmpEEG, ClustPars); % This is the actual clustering within subjects
    tmpEEG.msinfo.MSMaps(4).Maps = ThomasVersion'
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, tmpEEG, AllSubjects (1)); % Done, we just need to store this
% end

eeglab redraw

% save(fullfile(SavePathName,SaveFileName));
% % 
% % %% Let's see what a good number of clusters would be like
% % 
% % %pop_BootstrapMSNumber(ALLEEG,1:numel(ALLEEG), true, 50, 100);
% % 
% % %% Now we combine the microstate maps across subjects and sort the mean
% % 
% pluginpath = fileparts(which('eegplugin_Microstates.m'));                  % Getting eeglab path
% templatepath = fullfile(pluginpath,'Templates');
% 
% tmpEEG = pop_loadset('filename','Normative microstate template maps Neuroimage 2002.set','filepath',templatepath);
% [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
% 
%  NormativeTemplateIndex = CURRENTSET;
% % 
% % for Group = 1:nGroups
% %     % The mean of group X
% %     EEG = pop_CombMSTemplates(ALLEEG, GroupIndex{Group}, 0, 0, sprintf('GrandMean Group %i',Group));
% %     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
% %     [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
% %     GrandMeanIndex(Group) = CURRENTSET; % And keep track of it
% % end
% % 
% % % Now we want the grand-grand mean, based on the group means
% % if nGroups > 1
% %     EEG = pop_CombMSTemplates(ALLEEG, GrandMeanIndex, 1, 0, 'GrandGrandMean');
% %     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
% %     %?@We automatically sort the grandgrandmean based on a template from the literature
% %     GrandGrandMeanIndex = CURRENTSET; % and keep track of it
% % else
% %     GrandGrandMeanIndex = GrandMeanIndex(1);
% % end
% 
% GrandGrandMeanIndex 
% 
% [ALLEEG,EEG] = pop_SortMSTemplates(ALLEEG, 43, 1, 44);
% [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
% 
% % This should now be as good as possible, but we should look at it
% [ALLEEG,EEG] = pop_ShowIndMSMaps(EEG, 4, 1, ALLEEG); % Here, we go interactive to allow the user to put the classes in the canonical order
% [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % , we store it
% 
% 
% eeglab redraw
% save(fullfile(SavePathName,SaveFileName));
% 
% %% And we sort things out over means and subjects
% % First, the sequence of the two group means has be adjusted based on the
% % grand grand mean
% if nGroups > 1
%     ALLEEG = pop_SortMSTemplates(ALLEEG, GrandMeanIndex, 1, GrandGrandMeanIndex);
%     [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
% end
% 
% % Then, we sort the individuals based on their group means
% for Groups = 1:nGroups
%     ALLEEG = pop_SortMSTemplates(ALLEEG, GroupIndex{Group}, 0, GrandMeanIndex(Group)); % Group 1
%     [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
% 
% end
% 
% EEG = eeg_retrieve(ALLEEG,CURRENTSET);
% 
% eeglab redraw
% save(fullfile(SavePathName,SaveFileName));
% 
% %% eventually save things
% 
% for f = 1:numel(ALLEEG)
%     EEG = eeg_retrieve(ALLEEG,f);
%     fname = EEG.setname;
%     pop_saveset( EEG, 'filename',fname,'filepath',SavePath);
% end
% 
% %% Visualize some stuff to see if the fitting parameters appear reasonable
% % These are the paramters for the continuous fitting
% % FitPars = struct('nClasses',4,'lambda',1,'b',30,'PeakFit',false, 'BControl',true);
% 
% 
% 
% % Just a look at the first EEG
% EEG = eeg_retrieve(ALLEEG,1); 
% pop_ShowIndMSDyn([],EEG,0,FitPars);
% pop_ShowIndMSMaps(EEG,FitPars.nClasses);

%% Here comes the stats part

% Using the individual templates
%pop_QuantMSTemplates(ALLEEG, AllSubjects, 0, FitPars, []                   , fullfile(SavePath,'ResultsFromIndividualTemplates.csv'));

% And using the grand grand mean template
pop_QuantMSTemplates(ALLEEG, AllSubjects, 1, FitPars, 43, fullfile(SavePath,'ResultsFromGrandGrandMeanTemplate.csv'));

%% Eventually export the individual microstate maps to do statistics in Ragu

% nMaps = 4;
% 
% Grouping = nan(numel(AllSubjects),1);
% 
% for Groups = 1:nGroups
%     Grouping(GroupIndex{Group}) = Group;
% end
% 
% rd = SaveMSMapsForRagu(ALLEEG(AllSubjects),nMaps,Grouping);
% 
% save(fullfile(SavePath,'IndividualMaps.mat'),'rd');