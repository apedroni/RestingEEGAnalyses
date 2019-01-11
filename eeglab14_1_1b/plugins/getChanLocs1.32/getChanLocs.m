% getChanLocs(): Populate EEG.chanlocs using Wavefront .obj created by 3D-scanners (e.g. structure.io).
% FieldTrip toolbox functions are adapted to localize electrodes. See the original process at 
% http://www.fieldtriptoolbox.org/tutorial/electrode
%
% Usage:
%   >>EEG = getChanLocs(EEG, objPath, 'key1', value1, ..., 'keyN', valueN);
%
% Inputs:
%   EEG     - EEGLAB EEG structure
%   objPath - full filepath to folder containing [Model.obj, Model.jpg, Model.mtl] files
%
% Optional Inputs (value 0 or false to turn off):
%   'anonymizeFace'   - (Default = 1) Overwrite objPath .jpg file to replace 
%                        skintones with grey to anonymize subject's face
%   'chanLabels'      - (Default = {EEG.chanlocs(1,:).labels}') Label names
%                        for EEG channels (and misc sensors) to be localized.
%                        Default channel label list is extracted from EEG recording. 
%   'deleteTxtOutput' - (Default = 1) Delete text file containing electrode labels and 
%                        locations after importing to EEGLAB .set file.
%   'moveElecInwards' - (Default = 7.50) Move electrode locations towards (in mm)
%                        scalp to adjust for cap and electrode well thickness.
%                        Negative numbers result in an outward move, away from the scalp.
%   'saveName'        - (Default =  strcat(objPath, filesep, 'getChanLocs.txt')) Full filename 
%                        (including path) of output file containing electrode labels
%                        and locations. Imported into EEGLAB using readlocs(). Can be set to 
%                        automatically delete after import (see Optional Inputs: deleteTxtOutput)
%
% See also:
%   readlocs, ft_read_headshape, ft_electrodeplacement, ft_meshrealign
%   For FieldTrip (ft_*) functions, see http://www.fieldtriptoolbox.org
%   for the documentation and details.
%
% Author: 
%   Clement Lee, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego
%
% History: 
%   26 Jan 2018 v1.32 CL. try/catch on moveElecInwards for now so that process continues and locations are saved.
%   25 Jan 2018 v1.31 CL. Addressing issues with mex files when solid_angle.m (for moveElecInwards) fails. 
%   25 Jan 2018 v1.3  CL. Version check for string() when writing text file.
%   23 Jan 2018 v1.2  CL. Rename from GetChanLocs to getChanLocs. Default moveElecInwards 12.5-> 7.5.
%   19 Jan 2018	v1.1  CL. Switch fiducials to L/RHT (from L/RPA). Removed [] for channel / fiducial names. 
%   20 Dec 2017 v1.0  CL. Special thanks to M. Milham and L. Ai (@Child Mind Institute); N. Langer (@University of Zurich);
%					     and H. Tanaka (@Japan Advanced Institute of Science and Technology) for interest in testing.
%   4  Dec 2017 v0.1 Clement Lee. Created.
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function EEG = getChanLocs(EEG, objPath, varargin)
if nargin < 1
	help getChanLocs;
	return;
elseif mod(nargin,2) == 1
    help getChanLocs;
    error('Please check input format.')
elseif ~(any(size(dir([objPath filesep '*.obj']),1))||...
         any(size(dir([objPath filesep '*.jpg']),1))||...
         any(size(dir([objPath filesep '*.mtl']),1)))
    error('Please input path to folder containing [Model.obj, Model.jpg, Model.mtl] files.')
end

opts = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
if ~isfield(opts,'anonymizeFace')
    opts.anonymizeFace = 1; end
if ~isfield(opts,'chanLabels')
    opts.chanLabels = {EEG.chanlocs(1,:).labels}'; end
if ~isfield(opts,'deleteTxtOutput')
    opts.deleteTxtOutput = 1; end
if ~isfield(opts,'moveElecInwards')
    opts.moveElecInwards = 7.50; end
if ~isfield(opts,'renameAxis')
    opts.renameAxis = false; end 
if ~isfield(opts,'saveName')
    opts.saveName = strcat(objPath, filesep, 'getChanLocs.txt');
elseif isempty(regexp(opts.saveName, '.txt','once')) % endsWith(opts.saveName, '.txt')
    opts.saveName = strcat(opts.saveName,'.txt'); end
fprintf('Electrode location .txt file will be written to %s\n', opts.saveName)

if opts.anonymizeFace
fprintf('Anonymizing face by replacing skintones with grey...\n')
anonFace(strcat(objPath, filesep, 'Model.jpg'));
end

fprintf('Loading 3D model in mm scale...\n')
head_surface = ft_read_headshape(strcat(objPath, filesep, 'Model.obj'), 'unit','mm');

fprintf('Select (in order) Nasion, Left Helix/Tragus Intersection, and Right Helix/Tragus Intersection...\n')
cfg = [];
cfg.method = 'headshape';
cfg.channel = {'NAS','LHT','RHT'}';
fiducials = ft_electrodeplacement(cfg,head_surface);
close gcf

fprintf('Using fiducials to align to BTi coordinates...\n')
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'bti';
cfg.fiducial.nas    = fiducials.elecpos(1,:); %position of NAS
cfg.fiducial.lpa    = fiducials.elecpos(2,:); %position of LHT
cfg.fiducial.rpa    = fiducials.elecpos(3,:); %position of RHT
head_surface = ft_meshrealign(cfg,head_surface);

fprintf('Select electrode locations...\n')
cfg = [];
cfg.method = 'headshape';
cfg.channel = opts.chanLabels;
elec = ft_electrodeplacement(cfg,head_surface);
close gcf

if opts.moveElecInwards
	try
    cfg = [];
    cfg.method = 'moveinward';
    cfg.moveinward = opts.moveElecInwards;
    cfg.elec = elec;
    elec = ft_electroderealign(cfg);
	fprintf('Moving electrode locations from cap surface in towards scalp by %2.2f mm...\n', opts.moveElecInwards)
    catch
	warning('Failed to move electrodes inwards. Bug is being addressed. Proceeding to save locations as is...')
end

%% format and save ascii for import. delete file afterwards if requested
fprintf('Writing electrode locations to txt file...\n')
fileID = fopen(opts.saveName,'w');

v = ver('Matlab');
if str2double(v.Version)>=9.1
% 'labels','X','Y','Z'
fprintf(fileID,'%6s %9.4f %9.4f %9.4f\n', [string(elec.label) ; elec.elecpos(:,:)']);
else
    for ii = 1:length(elec.label)
    fprintf(fileID, '%6s %9.4f %9.4f %9.4f\n', elec.label{ii}, elec.elecpos(ii,:));
    end
end
fclose(fileID);

fprintf('Importing locations with readlocs()...\n')
EEG.chanlocs = readlocs(opts.saveName,'format',{'labels','X','Y','Z'});
if opts.deleteTxtOutput == 1
    fprintf('Deleting txt file with electrode locations...\n')
    delete(opts.saveName)
end
fprintf('Electrode Localization by 3D Object Finished!\n')
end
end

function anonFace(objJpg) %anonymize face by replacing skintones with grey
ogI = imread(objJpg);
% whiten all black points for processing
I = ogI; I(I==0) = 255;
% extract rgb as double
red  = double(I(:,:,1)); green = double(I(:,:,2)); blue = double(I(:,:,3));
% normalize by red and mask skintone values
greenR = green./red; blueR = blue./red;
gMask = (greenR>0.5 & greenR<0.85); uMask = (blueR>0.3 &  blueR<0.65);
mask = gMask & uMask;
ogI(mask(:,:,[1,1,1])) = 128;
% Overwrite Image
imwrite(ogI, objJpg)
end