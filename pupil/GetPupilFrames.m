function [croppedPupilFrames, thresh, eyeMask] = GetPupilFrames(dataBasePath, dataName, varargin)
%%% Function: Saves out configuration file with crop coordinates, threshold for segmenting the pupil, and eye mask
%%%           If specified, loads a previously generated configuration file and returns cropped video, thresh, eyeMask 
%%% Inputs:
%%% dataBasePath - full path to folder containing day of experiments i.e. media/data/
%%% dataName - specific experiment name, i.e. 20150101/vglut1_thinskull_13
%%% pupilTrackFromManta - if false, will track using Orca frame
%%% cropCoordFile - (optional) if provided, then will only set a new
%%%                 threshold value, but will maintain provided crop
%%%                 coordinates. Provide the dataName to use. 
%%% Output:
%%% croppedPupilFrames - the cropped pupil video. 
%%% thresh - pupil threshold value
%%% eyeMask - outline of the eyeball (i.e. maximum size pupil can achieve)
%%%
%%% If you need to just change the threshold after having already generated
%%% a config file here, use SetPupilThresh(dataName, newThresh).

global basePath

p = inputParser();
p.addParameter('templateDataName', '', @ischar); % Same format as dataName. If not empty, will load a previously generated pupil configuration file
% p.addParameter('setThresh', [], @isnumeric); % Keep the same crop coordinates, as indicated in the provided cropCoordFile, but change the threshold.
p.addParameter('pupilTrackFromOrca', false, @islogical); % Default to pupil tracking from manta
p.addParameter('isTwoColor', true, @islogical); % If using orca frames for pupil tracking, specify whether orca frames are two colors (i.e. 470 and 405). 
p.addParameter('fromBpod', false, @islogical); % if data was acquired with bpod
p.parse(varargin{:});


templateDataName = p.Results.templateDataName;
% setThresh = p.Results.setThresh;
pupilTrackFromOrca = p.Results.pupilTrackFromOrca;
isTwoColor = p.Results.isTwoColor;
fromBpod = p.Results.fromBpod;

loadMantaFrames = ~pupilTrackFromOrca;

if ~fromBpod
    mantaDir = fullfile(dataBasePath,[dataName,'_manta.mp4']);
    processedDataDir = fullfile(basePath, 'processedData', dataName);
else
%     dataFileName = ['Trial', num2str(dataName, '%05d')];
%     mantaDir = fullfile(dataBasePath, [dataFileName, '.mp4']);
    mantaDir = fullfile(dataBasePath, dataName);
    processedDataDir = strrep(dataBasePath, 'data', 'processedData');
    mkdir(processedDataDir);
end


if loadMantaFrames
    disp(['Loading manta dir ', mantaDir])
    [~, ~, mantaFrames] = LoadMP4(mantaDir);
    disp('Done loading manta')
end
if pupilTrackFromOrca
    imageDir = fullfile(dataBasePath,[dataName,'_orca']);
    if ~isTwoColor
        img470 = LoadImageStack( imageDir );
    else
        [img470, ~] = DualColorLoadVid(imageDir);
    end
    mantaFrames = img470;
end

generateNewConfiguration = isempty(templateDataName); % Recrop pupil etc. 

%%% Crop pupil
if generateNewConfiguration
    disp('Click top left and bottom right corner of crop region around eye')
    [croppedPupilFrames, x, y] = CropVid(mantaFrames);
    pupilInfo.x = x;
    pupilInfo.y = y;
else
    if fromBpod
        load(fullfile(templateDataName, 'pupilCropCoords.mat')); % Loads pupilInfo struct
    else
        load(fullfile(basePath, 'processedData', templateDataName, 'pupilCropCoords.mat')); % Loads pupilInfo struct
    end
     [croppedPupilFrames, x, y] = CropVid(mantaFrames, pupilInfo.x, pupilInfo.y);
end

%%% Now set pupil threshold. Click once in pupil, and once outside pupil.
if generateNewConfiguration 
     doAgain = true;
     while doAgain
        img = croppedPupilFrames(:,:,round(end/2));
        figure, imagesc(img);
%         [x, y] = ginput(1);
%         thresh = croppedPupilFrames(round(y), round(x), 1)
        res = input('Set threshold: ');
        thresh = res;
        img(img < thresh) = 0;
        pupilInfo.thresh = thresh;
        figure, imagesc(img);
        res = input('Do again (type Y then enter, else just enter)? ', 's');
        doAgain = isequal(upper(res),'Y');
     end
else
    thresh = pupilInfo.thresh;
end

%%% Set outline of eye. Draw polygon around edge of eye. 
disp('Now click a polygon around the white of the eye (excluding dark regions at edge that may be confused with pupil.')
disp('Double click in window when finished.')
if generateNewConfiguration
    eyeMask = roipoly(); %%% Outline the eye here.
    pupilInfo.eyeMask = eyeMask;
else
    eyeMask = pupilInfo.eyeMask;
end

saveFile = fullfile(processedDataDir, 'pupilCropCoords.mat');
save(saveFile, 'pupilInfo');

clear mantaFrames
disp('Done loading manta')