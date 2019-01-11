% eegplugin_Microstates() - EEGLAB plugin for microstate analyses
%
% Usage:
%   >> eegplugin_Microstates(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
%   Integrates the following subfunctions:
%   --------------------------------------
%   Computational tools:
%   <a href="matlab:helpwin pop_FindMSTemplates">pop_FindMSTemplates</a>  (identify microstate cluster maps)
%   <a href="matlab:helpwin pop_CombMSTemplates">pop_CombMSTemplates</a>  (average microstate cluster maps)
%   <a href="matlab:helpwin pop_SortMSTemplates">pop_SortMSTemplates</a>  (sorts microstate cluster maps based on a mean)
%   <a href="matlab:helpwin pop_QuantMSTemplates">pop_QuantMSTemplates</a> (quantifies the presence of microstates in EEG data)
%
%   Visualisations:
%   <a href="matlab:helpwin pop_ShowIndMSMaps">pop_ShowIndMSMaps</a>   (view and edit microstate cluster maps)
%   <a href="matlab:helpwin pop_ShowIndMSDyn">pop_ShowIndMSDyn</a>     (view microstate dynamics)
%
% A typical basic workflow looks like this:
% - Identify  microstate maps in each EEG (Tools menu)
% - Average microstate maps across the datasets (Tools menu)
% - Edit the mean microstate maps for the desired sequence (Plot menu)
% - Either
%     - Sort the individual microstate maps based on the edited average
%       (Tools menu)
%     - Quantify the microstate presence based on the sorted individual
%       microstate maps (Tools menu)
% - Or   
%     - Quantify the microstate presence based on the averaged
%       microstate maps (Tools menu)
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
function vers = eegplugin_Microstates (fig, try_strings, catch_strings)

    vers = 'Microstates0.5.0';

    global MSTEMPLATE;
    
    comFindMSTemplates     = [try_strings.no_check '[       EEG LASTCOM] = pop_FindMSTemplates(EEG);'                       catch_strings.store_and_hist]; % ok
    comCombineMSTemplates  = [try_strings.no_check '[       EEG LASTCOM] = pop_CombMSTemplates(ALLEEG,CURRENTSET,false);'   catch_strings.new_and_hist];
    comCombineMSMeans      = [try_strings.no_check '[       EEG LASTCOM] = pop_CombMSTemplates(ALLEEG,CURRENTSET,true );'   catch_strings.new_and_hist];
    comSortMSTemplates     = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_SortMSTemplates(ALLEEG,CURRENTSET,false);'   catch_strings.store_and_hist];
    comSortMSMeans         = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_SortMSTemplates(ALLEEG,CURRENTSET,true) ;'   catch_strings.store_and_hist];
    comSortMSTemplatesT    = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_SortMSTemplates(ALLEEG,CURRENTSET,false,MSTEMPLATE);'   catch_strings.store_and_hist];
    comSortMSMeansT        = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_SortMSTemplates(ALLEEG,CURRENTSET,true, MSTEMPLATE) ;'   catch_strings.store_and_hist];
    
    comGetIndMSDynamics    = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_GetMSDynamics(ALLEEG,EEG,false);'            catch_strings.new_and_hist];
    comGetMeanMSDynamics   = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_GetMSDynamics(ALLEEG,EEG,true);'             catch_strings.new_and_hist];
    comQuantMSTemplatesS   = [try_strings.no_check '[           LASTCOM] = pop_QuantMSTemplates(ALLEEG,CURRENTSET,false);'  catch_strings.add_to_hist];
    comQuantMSTemplatesM   = [try_strings.no_check '[           LASTCOM] = pop_QuantMSTemplates(ALLEEG,CURRENTSET, true);'  catch_strings.add_to_hist];
    comRaguMSTemplates     = [try_strings.no_check '[           LASTCOM] = pop_RaguMSTemplates(ALLEEG,CURRENTSET       );'  catch_strings.add_to_hist];
    %    comBootStrapMSNumber   = [try_strings.no_check '[           LASTCOM] = pop_BootstrapMSNumber(ALLEEG,CURRENTSET);'       catch_strings.add_to_hist];
    comSilhouetteMSNumber  = [try_strings.no_check '[           LASTCOM] = pop_MS_Silhouette(ALLEEG,CURRENTSET);'           catch_strings.add_to_hist];
    
    comShowIndMSMaps       = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_ShowIndMSMaps(EEG,[],false ,ALLEEG);'        catch_strings.add_to_hist];
    comEditIndMSMaps       = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_ShowIndMSMaps(EEG,[],true  ,ALLEEG);'        catch_strings.store_and_hist];
    comShowIndMSDyn        = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_ShowIndMSDyn(ALLEEG,EEG,false);'             catch_strings.store_and_hist];
    comShowIndMSDynM       = [try_strings.no_check '[ALLEEG EEG LASTCOM] = pop_ShowIndMSDyn(ALLEEG,EEG,true);'              catch_strings.store_and_hist];

    toolsmenu = findobj(fig, 'tag', 'tools');
    toolssubmenu = uimenu( toolsmenu, 'label', 'Microstates','userdata','study:on');

    plotmenu = findobj(fig, 'tag', 'plot');
    uimenu( plotmenu, 'label', 'Plot microstate maps'                       ,'CallBack',comShowIndMSMaps,'Separator','on');
    uimenu( plotmenu, 'label', 'Edit microstate maps'                       ,'CallBack',comEditIndMSMaps,'Separator','off');
    uimenu( plotmenu, 'label', 'Plot microstate dynamics'                   ,'CallBack',comShowIndMSDyn,'Separator','on');
    uimenu( plotmenu, 'label', 'Plot microstate dynamics (Mean template)'   ,'CallBack',comShowIndMSDynM);

    % create menus if necessary
    % -------------------------
    uimenu( toolssubmenu, 'Label', 'Identify  microstate  maps'                           , 'CallBack', comFindMSTemplates,    'position', 1);

    uimenu( toolssubmenu, 'Label', 'Compute mean microstate maps across individual'        , 'CallBack', comCombineMSTemplates, 'position', 2, 'userdata', 'study:on','Separator','on');
    uimenu( toolssubmenu, 'Label', 'Compute grand mean microstate maps across means'       , 'CallBack', comCombineMSMeans    , 'position', 3, 'userdata', 'study:on');
%    uimenu( toolssubmenu, 'Label', 'Crossvalidate optimal number of classes'              , 'CallBack', comBootStrapMSNumber , 'position', 4, 'userdata', 'study:on');
    uimenu( toolssubmenu, 'Label', 'Use silhouette plots to choose number of classes'      , 'CallBack', comSilhouetteMSNumber , 'position', 4, 'userdata', 'study:on');    

    uimenu( toolssubmenu, 'Label', 'Sort individual microstate maps according to mean'     , 'CallBack', comSortMSTemplates, 'position', 5, 'userdata', 'study:on','Separator','on');
    uimenu(toolssubmenu, 'label',  'Sort individual microstate maps according to template' , 'CallBack', comSortMSTemplatesT,'position', 6, 'userdata', 'study:on');

    uimenu( toolssubmenu, 'Label', 'Sort mean microstate maps accoring to a grand mean'    , 'CallBack', comSortMSMeans    , 'position', 7, 'userdata', 'study:on','Separator','on');    
    uimenu(toolssubmenu, 'label',  'Sort mean microstate maps according to template'       , 'CallBack', comSortMSMeansT,    'position',8, 'userdata', 'study:on');

    uimenu( toolssubmenu, 'Label', 'Obtain microstate dynamics (own template maps)'  , 'CallBack', comGetIndMSDynamics    ,   'position',  9, 'Separator','on');
    uimenu( toolssubmenu, 'Label', 'Obtain microstate dynamics (mean template maps)' , 'CallBack', comGetMeanMSDynamics   ,   'position',  10);

    uimenu( toolssubmenu, 'Label', 'Quantify microstates in dataset (own template maps)'  , 'CallBack', comQuantMSTemplatesS,   'position', 11, 'Separator','on');
    uimenu( toolssubmenu, 'Label', 'Quantify microstates in dataset (mean template maps)' , 'CallBack', comQuantMSTemplatesM,   'position', 12);
    if numel(which('Ragu')) > 0
        uimenu( toolssubmenu, 'Label', 'Test for topographic effects in microstate topographies (Ragu)' , 'CallBack', comRaguMSTemplates,   'position', 13,'Separator','on');
    end
    pluginpath = fileparts(which('eegplugin_Microstates.m'));                  % Get eeglab path
    templatepath = fullfile(pluginpath,'Templates');

    Templates = dir(fullfile(templatepath,'*.set'));
       
    for t = 1: numel(Templates)
        MSTemplate(t) = pop_loadset('filename',Templates(t).name,'filepath',templatepath);
    end
    
    MSTEMPLATE = MSTemplate;

end