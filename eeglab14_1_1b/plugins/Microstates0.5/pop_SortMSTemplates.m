%pop_SortMSTemplates() Reorder microstate maps based on a mean template
%
% Usage: >> [AllEEG,EEGout,com] = pop_SortMSTemplates(AllEEG,SetToSort, DoMeans, TemplateSet, IgnorePolarity)
%
% EEG lab specific:
%
%   "AllEEG" 
%   -> AllEEG structure with all the EEGs that may be analysed
%   
%   EEGout
%   -> The refreshed set of CURRENTSET EEGs
%
%   "SetToSort" 
%   -> Index of selected EEGs. If more than one EEG is selected, the analysis
%      will be limited to those, if not, the user is asked.
%
%   "TemplateSet"
%   -> Index of the AllEEG element with the dataset used as a template for sorting
%   If TemplateSet is -1, the user will be prompted for a normative
%   template. If you want to use a normative template in a script, load the
%   dataset with the normative template first and than make this the
%   TemplateSet that you use for sorting. If the template set is an EEG
%   structure, or an array of EEG structures, these will be used.
%
%   "IgnorePolarity"
%   -> Ignore the polarity of the maps to be sorted   
%
% Output:
%
%   "AllEEG" 
%   -> AllEEG structure with all the updated EEGs
%
%   "com"
%   -> Command necessary to replicate the computation
%
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
%
function [AllEEG,EEGout,com] = pop_SortMSTemplates(AllEEG,SetToSort, DoMeans, TemplateSet, IgnorePolarity)

    com = '';
    
%    if numel(EEG) > 1
%        errordlg2('pop_findMSTemplates() currently supports only a single EEG as input');
%        return;
%    end
    SelectedSets = SetToSort;
    if nargin < 3  DoMeans = false;            end %#ok<*SEPEX>
    if nargin < 5  IgnorePolarity = true;      end 
    
    nonempty = find(cellfun(@(x) isfield(x,'msinfo'), num2cell(AllEEG)));
    HasChildren = cellfun(@(x) isfield(x,'children'), {AllEEG.msinfo});
    nonemptyInd  = nonempty(~HasChildren);
    nonemptyMean = nonempty(HasChildren);
    
%    if numel(nonemptyMean) < 1
%        errordlg2('No mean templates found','Sort microstate classes');
%        return;
%    end
    
    if nargin < 4
        TemplateSet = nonemptyMean(1); 
    end
    
    if isstruct(TemplateSet)
        TemplateNames = {TemplateSet.setname};
        MeanIndex = 1;
    else
        TemplateNames = {AllEEG(nonemptyMean).setname};
        MeanIndex = find(nonemptyMean == TemplateSet,1);
    end
           
  
    if numel(SetToSort) < 2
        if DoMeans == true 
            AvailableSets = {AllEEG(nonemptyMean).setname};
        else
            AvailableSets = {AllEEG(nonemptyInd).setname};
        end
                  
        res = inputgui( 'geometry', {1 1 1 1 1 1}, 'geomvert', [1 1 4 1 1 1], 'uilist', { ...
            { 'Style', 'text', 'string', 'Choose sets for sorting'} ...
            { 'Style', 'text', 'string', 'Use ctrlshift for multiple selection'} ...
            { 'Style', 'listbox', 'string', AvailableSets, 'tag','SelectSets' ,'Min', 0, 'Max',2} ...
            { 'Style', 'text', 'string', 'Name of mean', 'fontweight', 'bold'  } ...
            { 'Style', 'popupmenu', 'string', TemplateNames,'tag','MeanName','Value',MeanIndex} ...
            { 'Style', 'checkbox', 'string', 'No polarity','tag','Ignore_Polarity' ,'Value', IgnorePolarity }  ...
            });
     
        if isempty(res); return; end
    
        MeanIndex = res{2};
        
        IgnorePolarity = res{3};
        if DoMeans == true
            SelectedSet = nonemptyMean(res{1});
        else
            SelectedSet = nonemptyInd(res{1});
        end
    else
        if nargin < 4
            res = inputgui( 'geometry', {1 1 1}, 'geomvert', [1 1 1], 'uilist', { ...
                { 'Style', 'text', 'string', 'Name of mean', 'fontweight', 'bold'  } ...
                { 'Style', 'popupmenu', 'string', TemplateNames,'tag','MeanName' ,'Value',MeanIndex } ...
                { 'Style', 'checkbox', 'string', 'No polarity','tag','Ignore_Polarity' ,'Value', IgnorePolarity }  ...
                });
        
            if isempty(res); return; end
            MeanIndex = res{1};
            IgnorePolarity = res{2};
        end    
        SelectedSet = SetToSort;
    end

    if isstruct(TemplateSet)
        ChosenTemplate = TemplateSet(MeanIndex);
    else
        ChosenTemplate = AllEEG(1,nonemptyMean(MeanIndex));
    end
    
    if numel(SelectedSet) < 1
        errordlg2('You must select at least one set of microstate maps','Sort microstate classes');
        return
    end

    MinClasses     = ChosenTemplate.msinfo.ClustPar.MinClasses;
    MaxClasses     = ChosenTemplate.msinfo.ClustPar.MaxClasses;
%    IgnorePolarity = ChosenTemplate.msinfo.ClustPar.IgnorePolarity;
%    GFPPeaks       = ChosenTemplate.msinfo.ClustPar.GFPPeaks;
    
    for index = 1:length(SelectedSet)
        sIndex = SelectedSet(index);
        if ~isfield(AllEEG(sIndex),'msinfo')
            errordlg2(sprintf('Microstate info not found in datset %s',AllEEG(sIndex).setname),'Sort microstate classes'); 
            return;
        end
    
        if  MinClasses     < AllEEG(sIndex).msinfo.ClustPar.MinClasses || ...
            MaxClasses     > AllEEG(sIndex).msinfo.ClustPar.MaxClasses
%            IgnorePolarity ~= AllEEG(sIndex).msinfo.ClustPar.IgnorePolarity || ...
%            GFPPeaks       ~= AllEEG(sIndex).msinfo.ClustPar.GFPPeaks;
                errordlg2('Microstate parameters differ between datasets','Sort microstate classes');
                return;
        end
    end

    for n = MinClasses:MaxClasses
        MapsToSort = nan(numel(SelectedSet),n,numel(ChosenTemplate.chanlocs));
        % Here we go to the common set of channels
        for index = 1:length(SelectedSet)
            sIndex = SelectedSet(index);
            LocalToGlobal = MakeResampleMatrices(AllEEG(sIndex).chanlocs,ChosenTemplate.chanlocs);
           MapsToSort(index,:,:) = AllEEG(sIndex).msinfo.MSMaps(n).Maps * LocalToGlobal';
        end
        % We sort out the stuff
        [~,SortOrder, Communality, polarity] = ArrangeMapsBasedOnMean(MapsToSort,ChosenTemplate.msinfo.MSMaps(n).Maps,~IgnorePolarity);
        %   And we go back to the original channel order
         for index = 1:length(SelectedSet)
            sIndex = SelectedSet(index);
            AllEEG(sIndex).msinfo.MSMaps(n).Maps = AllEEG(sIndex).msinfo.MSMaps(n).Maps(SortOrder(index,:),:);
            AllEEG(sIndex).msinfo.MSMaps(n).Maps = AllEEG(sIndex).msinfo.MSMaps(n).Maps .* repmat(polarity(index,:)',1,numel(AllEEG(sIndex).chanlocs));
            AllEEG(sIndex).msinfo.MSMaps(n).ColorMap = ChosenTemplate.msinfo.MSMaps(n).ColorMap;
            AllEEG(sIndex).msinfo.MSMaps(n).SortMode = 'template based';
            if strcmp(AllEEG(sIndex).msinfo.MSMaps(n).SortedBy,'none')
                AllEEG(sIndex).msinfo.MSMaps(n).SortedBy = ChosenTemplate.setname;
            else
                AllEEG(sIndex).msinfo.MSMaps(n).SortedBy = [ChosenTemplate.msinfo.MSMaps(n).SortedBy '->' ChosenTemplate.setname];
            end
            AllEEG(sIndex).msinfo.MSMaps(n).Communality = Communality(index,:);
            AllEEG(sIndex).saved = 'no';
         end
    end

    txt = sprintf('%i ',SelectedSet);
    txt(end) = [];
    if isstruct(TemplateSet)
        txt2 = inputname(4);
        if isempty(txt2)
            txt2 = '[]';
        end
    else
        txt2 = sprintf('%i ',TemplateSet);
        txt2(end) = [];
    end
    
    EEGout = AllEEG(SelectedSets);
    
    com = sprintf('[%s %s com] = pop_SortMSTemplates(%s, [%s], %i, %s, %i);', inputname(1),inputname(2),inputname(1), txt, DoMeans, txt2,IgnorePolarity);
end
