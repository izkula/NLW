function [vids, vidnames, viddirs] = Load2pImageTrials(imageDir, dataName, varargin)
% Load all bpod trials from a session, saved in the directory specified
% by imageDir i.e. 
%   '/media/Data/data/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2'
%
% Returns: A cell array, where each entry is the video for each trial
%          vidnames: cell array where each entry is path to first tif file
%          of each trial
%          viddirs: cell array where each entry is the path to the folder
%          containing the tif files for a trial

    p = inputParser();
    p.addParameter('subSeq', [], @isnumeric);
    p.addParameter('doLoad', true, @islogical)
    p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
    p.addParameter('zplane', [], @isnumeric);
    
    p.parse(varargin{:});
    subSeq = p.Results.subSeq;
    doLoad = p.Results.doLoad;
    isNLW = p.Results.isNLW;
    zplane = p.Results.zplane;

    trials = dir(imageDir);
    numTrials = numel(trials);
    if isempty(subSeq)
        subSeq = 1:(numTrials-3); %%% Exclude the final session and the '.' and '..' files
    end
    subSeq = subSeq(subSeq <= (numTrials - 3));

    trialNames = {};
    for i=1:numTrials
        if ~strcmp(trials(i).name, '.') && ~strcmp(trials(i).name, '..')
            trialIdx = strsplit(trials(i).name,'Trial');
            if numel(trialIdx) > 1
                trialIdx = str2num(trialIdx{2});
                trialNames{trialIdx} = trials(i).name;
            end
        end
    end

    vids = cell(numel(subSeq), 1);
    viddirs = {};
    vidnames = {};
    for kk=subSeq
%         kk
        if isNLW
            if isempty(zplane)
                viddirs{kk} = fullfile(imageDir, trialNames{kk}, [dataName, '_', trialNames{kk}] );
            else
                viddirs{kk} = fullfile(imageDir, trialNames{kk}, [dataName, '_', trialNames{kk}], ['z', num2str(zplane)] );
            end
        else
            viddirs{kk} = fullfile(imageDir, trialNames{kk}, [dataName, '-', trialNames{kk}(end-2:end)] );
        end
        
        [vids{kk}, vidnames{kk}] = LoadImageStack( viddirs{kk}, [], doLoad, 'isNLW', isNLW, 'zplane', zplane);
    end
end
