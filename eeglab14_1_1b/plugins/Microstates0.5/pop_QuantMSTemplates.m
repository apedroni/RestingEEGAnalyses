%pop_QuantMSTemplates() quantifies the presence of microstates in EEG data
%
% Usage:
%   >> com = pop_QuantMSTemplates(AllEEG, CURRENTSET, UseMeanTmpl, FitParameters, MeanSet, FileName)
%
% EEG lab specific:
%
%   "AllEEG" 
%   -> AllEEG structure with all the EEGs that may be analysed
%
%   "CURRENTSET" 
%   -> Index of selected EEGs. If more than one EEG is selected, the analysis
%      will be limited to those, if not, the user is asked
%
% Graphical interface / input parameters
%
%   UseMeanTmpl
%   -> True if a mean cluster center is to be used to quantify the EEG
%   data, false if the template from the data itself is to be used
%
%   FitParameters 
%   -> A struct with the following parameters:
%      - nClasses: The number of classes to fit
%      - PeakFit : Whether to fit only the GFP peaks and interpolate in
%        between (true), or fit to the entire data (false)
%      - b       : Window size for label smoothing (0 for none)
%      - lambda  : Penalty function for non-smoothness
%
%   "Name of Mean" (only for the GUI)
%   -> EEG dataset containing the mean clusters to be used if UseMeanTmpl
%   is true, else not relevant
%
%   Meanset
%   -> Index of the AllEEG dataset containing the mean clusters to be used if UseMeanTmpl
%   is true, else not relevant
%   Filename
%   -> Name of the file to store the output. You can store CSV files of
%   Matlab files
%
% Output:
%
%   "com"
%   -> Command necessary to replicate the computation
%              %
% Author: Thomas Koenig, University of Bern, Switzerland, 2016
%
% Copyright (C) 2016 Thomas Koenig, University of Bern, Switzerland, 2016
% thomas.koenig@puk.unibe.ch
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function com = pop_QuantMSTemplates(AllEEG, CURRENTSET, UseMeanTmpl, FitParameters, MeanSet, FileName)

    if nargin < 2,  CURRENTSET    = [];     end   
    if nargin < 3,  UseMeanTmpl   = false;  end
    if nargin < 4,  FitParameters = [];     end 
    if nargin < 5,  MeanSet       = [];     end 

    com = '';
    if nargin < 3
        ButtonName = questdlg('What type of templates do  you want to use?', ...
                         'Microstate statistics', ...
                         'Sorted individual maps', 'Averaged maps', 'Sorted individual maps');
        switch ButtonName,
            case 'Individual maps',
                UseMeanTmpl = false;
            case 'Averaged maps',
                UseMeanTmpl = true;
        end % switch
    end 
    
    nonempty = find(cellfun(@(x) isfield(x,'msinfo'), num2cell(AllEEG)));
    HasChildren = cellfun(@(x) isfield(x,'children'), {AllEEG.msinfo});
    nonemptyInd  = nonempty(~HasChildren);
    nonemptyMean = nonempty(HasChildren);
    
    AvailableMeans = {AllEEG(nonemptyMean).setname};
    AvailableSets  = {AllEEG(nonemptyInd).setname};
    
    if numel(CURRENTSET) > 1
        SelectedSet = CURRENTSET;
    
        if UseMeanTmpl == true && isempty(MeanSet) 
            disp('UseMeanTemplate')
            res = inputgui( 'geometry', {1 1}, 'geomvert', [1 1], 'uilist', { ...
                { 'Style', 'text', 'string', 'Name of mean', 'fontweight', 'bold'  } ...
                 { 'Style', 'listbox', 'string', AvailableMeans, 'tag','SelectSets'}});
     
            if isempty(res)
                return
            end
            MeanSet = nonemptyMean(res{1});
        end
    else
        if UseMeanTmpl == true && isempty(MeanSet)
            res = inputgui( 'geometry', {1 1 1 1 1 1}, 'geomvert', [1 1 4 1 1 4], 'uilist', { ...
                { 'Style', 'text', 'string', 'EEGs to analyze', 'fontweight', 'bold' } ...    
                { 'Style', 'text', 'string', 'Use ctrlshift for multiple selection'} ...
                { 'Style', 'listbox', 'string', AvailableSets, 'tag','SelectSets' ,'Min', 0, 'Max',2} ...
                { 'Style', 'text', 'string', ''} ...
                { 'Style', 'text', 'string', 'Name of mean', 'fontweight', 'bold'  } ...
                { 'Style', 'listbox', 'string', AvailableMeans, 'tag','SelectSets'}});
     
            if isempty(res)
                return
            end
            SelectedSet = nonemptyInd(res{1});
            MeanSet     = nonemptyMean(res{2});
        else
           res = inputgui( 'geometry', {1 1 1}, 'geomvert', [1 1 4], 'uilist', { ...
                { 'Style', 'text', 'string', 'EEGs to analyze', 'fontweight', 'bold' } ...    
                { 'Style', 'text', 'string', 'Use ctrlshift for multiple selection'} ...
                { 'Style', 'listbox', 'string', AvailableSets, 'tag','SelectSets' ,'Min', 0, 'Max',2} ...
                });
     
            if isempty(res)
                return
            end
            SelectedSet = nonemptyInd(res{1});
        end
    end
    
    AllSets = [SelectedSet MeanSet];
    MinClasses = max(cellfun(@(x) AllEEG(x).msinfo.ClustPar.MinClasses,num2cell(AllSets)));
    MaxClasses = min(cellfun(@(x) AllEEG(x).msinfo.ClustPar.MaxClasses,num2cell(AllSets)));
    
    if UseMeanTmpl == false
        if isfield(AllEEG(SelectedSet(1)).msinfo,'FitPar');     par = AllEEG(SelectedSet(1)).msinfo.FitPar;
        else par = [];
        end
    else
        if isfield(AllEEG(MeanSet).msinfo,'FitPar');            par = AllEEG(MeanSet).msinfo.FitPar;
        else par = [];
        end
    end
    
    [par,paramsComplete] = UpdateFitParameters(FitParameters,par,{'nClasses','lambda','PeakFit','b', 'BControl'});
 
    if ~paramsComplete
        par = SetFittingParameters(MinClasses:MaxClasses,par);
    end
    
%    MSStats = table();
    
    h = waitbar(0);
    set(h,'Name','Quantifying microstates, please wait...');
    set(findall(h,'type','text'),'Interpreter','none');

%    MSStats(numel(SelectedSet)).DataSet = '';
    
    for s = 1:numel(SelectedSet)
        sIdx = SelectedSet(s);
        waitbar((s-1) / numel(SelectedSet),h,sprintf('Working on %s',AllEEG(sIdx).setname),'Interpreter','none');
        DataInfo.subject   = AllEEG(sIdx).subject;
        DataInfo.group     = AllEEG(sIdx).group;
        DataInfo.condition = AllEEG(sIdx).condition;
        DataInfo.setname   = AllEEG(sIdx).setname;
        
        if UseMeanTmpl == false
            Maps = AllEEG(sIdx).msinfo.MSMaps(par.nClasses).Maps;
            AllEEG(sIdx).msinfo.FitPar = par;
            [MSClass,gfp,ExpVar] = AssignMStates(AllEEG(sIdx),Maps,par,AllEEG(sIdx).msinfo.ClustPar.IgnorePolarity);
            if ~isempty(MSClass)
 %              MSStats = [MSStats; QuantifyMSDynamics(MSClass,AllEEG(sIdx).msinfo,AllEEG(sIdx).srate, DataInfo, '<<own>>')];
                MSStats(s) = QuantifyMSDynamics(MSClass,gfp,AllEEG(sIdx).msinfo,AllEEG(sIdx).srate, DataInfo, '<<own>>',ExpVar);
            end
        else
            Maps = AllEEG(MeanSet).msinfo.MSMaps(par.nClasses).Maps;
            AllEEG(sIdx).msinfo.FitPar = par;
            LocalToGlobal = MakeResampleMatrices(AllEEG(sIdx).chanlocs,AllEEG(MeanSet).chanlocs);
            [MSClass,gfp,ExpVar] = AssignMStates(AllEEG(sIdx),Maps,par, AllEEG(MeanSet).msinfo.ClustPar.IgnorePolarity, LocalToGlobal);
            if ~isempty(MSClass)
%                MSStats = [MSStats; QuantifyMSDynamics(MSClass,AllEEG(sIdx).msinfo,AllEEG(sIdx).srate, DataInfo, AllEEG(MeanSet).setname)]; 
                MSStats(s) = QuantifyMSDynamics(MSClass,gfp,AllEEG(sIdx).msinfo,AllEEG(sIdx).srate, DataInfo, AllEEG(MeanSet).setname, ExpVar);
            end
        end
    end
    close(h);
    idx = 1;
    if nargin < 6
        [FName,PName,idx] = uiputfile({'*.csv','Comma separated file';'*.csv','Semicolon separated file';'*.txt','Tab delimited file';'*.mat','Matlab Table'; '*.csv','Text file for R'},'Save microstate statistics');
        FileName = fullfile(PName,FName);
    else
        if ~isempty(strfind(FileName,'.mat'))
            idx = 2;
        end
    end

%       writetable(MSStats,FileName);

    switch idx
        case 1
            SaveStructToTable(MSStats,FileName,',');
        case 2
            SaveStructToTable(MSStats,FileName,';');
        case 3
            SaveStructToTable(MSStats,FileName,sprintf('\t'));

        case 4
            save(FileName,'MSStats');
        case 5
            error('This is still under construction');
%            SaveStructToR(MSStats,FileName);
    end
    
    txt = sprintf('%i ',SelectedSet);
    txt(end) = [];

    com = sprintf('com = pop_QuantMSTemplates(%s, [%s], %i, %s, %i, ''%s'');', inputname(1), txt, UseMeanTmpl, struct2String(par), MeanSet, FileName);
end
