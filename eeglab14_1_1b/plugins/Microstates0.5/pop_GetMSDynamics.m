%pop_pop_GetMSDynamics() Show microstate dynamics over time
%
% Usage:
%   >> [AllEEG, TheEEG, com] = pop_GetMSDynamics(AllEEG,TheEEG,UseMean,FitPar, MeanSet)
%
% EEG lab specific:
%
%   "AllEEG" 
%   -> AllEEG structure with all the EEGs that may be analysed
%
%   "TheEEG" 
%   -> EEG structure with the EEG to search for templates
%
%   UseMean
%   -> True if a mean cluster center is to be used to quantify the EEG
%   data, false (default) if the template from the data itself is to be used
%
%   "Number of Classes" / FitPar.nClasses 
%   -> Number of clusters to quantify
%
%   "Name of Mean" (GUI)
%   -> EEG dataset containing the mean clusters to be used if UseMean is
%   true, else not relevant
%
%   Meanset (parameter)
%   -> Index of the AllEEG dataset containing the mean clusters to be used 
%      if UseMean is true, else not relevant
%
% Output:
%
%   "AllEEG" 
%   -> AllEEG structure with all the EEG (fitting parameters may be updated)
%
%   "TheEEG" 
%   -> EEG structure with the EEG (fitting parameters may be updated)
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

function [AllEEG, TheEEG, com] = pop_GetMSDynamics(AllEEG,TheEEG,UseMean,FitPar, MeanSet)

    com = '';
    
    
    if numel(TheEEG) > 1
        errordlg2('pop_findMSTemplates() currently supports only a single EEG as input');
        return;
    end
    
    if nargin < 3,     UseMean =  false;    end
    if nargin < 4,     FitPar  = [];        end 
    if nargin < 5,     MeanSet = [];        end 
    
    if UseMean == false && ~isfield(TheEEG,'msinfo')
        errordlg2('The data does not contain microstate maps','Show microstate dynamics');
        return;
    end

    if UseMean == true && isempty(MeanSet)
        nonempty = find(cellfun(@(x) isfield(x,'msinfo'), num2cell(AllEEG)));
        HasChildren = cellfun(@(x) isfield(x,'children'), {AllEEG.msinfo});
        nonemptyMean = nonempty(HasChildren);
    
        AvailableMeans = {AllEEG(nonemptyMean).setname};
        res = inputgui( 'geometry', {1 1}, 'geomvert', [1 4], 'uilist', { ...
            { 'Style', 'text', 'string', 'Name of mean', 'fontweight', 'bold'  } ...
            { 'Style', 'listbox', 'string', AvailableMeans, 'tag','SelectSets'} ...
            });
     
        if isempty(res)
            return
        end
        MeanSet = nonemptyMean(res{1});
    end
    
    if UseMean == false
        if isfield(TheEEG.msinfo,'FitPar');                params = TheEEG.msinfo.FitPar;
        else params = []; %#ok<SEPEX>
        end
    else
        if isfield(AllEEG(MeanSet).msinfo,'FitPar');     params = AllEEG(MeanSet).msinfo.FitPar;
       else params = []; %#ok<SEPEX>
        end
    end
    [FitPar,paramsComplete] = UpdateFitParameters(FitPar,params,{'nClasses','lambda','PeakFit','b','BControl','Normalize','Rectify'});

    if nargin < 4 || paramsComplete == false
        if UseMean == false     
            FitPar = SetFittingParameters(TheEEG.msinfo.ClustPar.MinClasses:TheEEG.msinfo.ClustPar.MaxClasses,FitPar, true);
        else
            FitPar = SetFittingParameters(AllEEG(MeanSet).msinfo.ClustPar.MinClasses:AllEEG(MeanSet).msinfo.ClustPar.MaxClasses,FitPar, true);
        end
    end
    if isempty(FitPar.nClasses);   return; end
    
    if UseMean == false
        Maps   = TheEEG.msinfo.MSMaps(FitPar.nClasses).Maps;
    else
        Maps   = AllEEG(MeanSet).msinfo.MSMaps(FitPar.nClasses).Maps;
        AllEEG(MeanSet).msinfo.FitPar = FitPar;
    end
    TheEEG.msinfo.FitPar = FitPar;

    if UseMean == true
        LocalToGlobal = MakeResampleMatrices(TheEEG.chanlocs,AllEEG(MeanSet).chanlocs);
        [MSClass,gfp,fit] = AssignMStates(TheEEG,Maps,FitPar,AllEEG(MeanSet).msinfo.ClustPar.IgnorePolarity,LocalToGlobal);
    else
        [MSClass,gfp,fit] = AssignMStates(TheEEG,Maps,FitPar,TheEEG.msinfo.ClustPar.IgnorePolarity);
    end
    
    if isempty(MSClass)
        return;
    end
    
    nChannels = TheEEG.nbchan;
    TheEEG.nbchan = FitPar.nClasses;
    newEEG = zeros(TheEEG.nbchan,TheEEG.pnts,TheEEG.trials);
 
    for i = 1:FitPar.nClasses
        chanlocs(i) = struct('labels',sprintf('MS-Class-%i',i),'type','','theta', nan,'radius',nan,'X',nan,'Y',nan,'Z',nan,'sph_theta',nan,'sph_phi',nan,'sph_radius',nan,'urchan',nan,'ref','none');
        for s = 1:TheEEG.trials
            assigned = MSClass(:,s) == i;
            if FitPar.Normalize == true
                newEEG(i,:,s) = assigned;
            else
                newEEG(i,assigned,s) = Maps(i,:) * TheEEG.data(:,assigned,s);
                if FitPar.Rectify == true
                    newEEG(i,assigned,s) = abs(newEEG(i,assigned,s));
                end
            end
        end
    end
    TheEEG.chanlocs = chanlocs;
    TheEEG.data = newEEG;

    if UseMean < 2
        com = sprintf('[AllEEG EEG com] = pop_GetMSDynamics(%s, %s, 0, %s);'    , inputname(1), inputname(2),struct2String(FitPar));
    else
        com = sprintf('[AllEEG EEG com] = pop_GetMSDynamics(%s, %s, 1, %s, %i);', inputname(1), inputname(2), struct2String(FitPar), MeanSet);
    end
end

