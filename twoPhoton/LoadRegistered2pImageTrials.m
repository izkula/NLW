function [vids, vidnames, viddirs] = LoadRegistered2pImageTrials(imageDir, dataName, varargin)
% Load all registered bpod trials from a session, saved in the processed data
% directory corresponding to imageDir, where imageDir is i.e.
%   '/media/Data/data/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2'
%
% Returns: A cell array, where each entry is the video for each trial
%          vidnames: cell array where each entry is path to first tif file
%          of each trial
%          viddirs: cell array where each entry is the path to the folder
%          containing the tif files for a trial

    global basePath bpodImagePath

    p = inputParser();
    p.addParameter('subSeq', [], @isnumeric);
    p.addParameter('doLoad', true, @islogical);
    p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
    p.addParameter('zplane', [], @isnumeric);

    
    p.parse(varargin{:});
    subSeq = p.Results.subSeq;
    doLoad = p.Results.doLoad;
    isNLW = p.Results.isNLW;
    zplane = p.Results.zplane;
    
    
    imageDir = GetProcessedDataDir2p(imageDir, basePath, bpodImagePath);

    
    trials = dir(imageDir);
    numTrials = numel(trials);
    if isempty(subSeq)
        subSeq = 1:(numTrials-3);  %%% Exclude the final session and the '.' and '..' files
    end
    subSeq = subSeq(subSeq <= (numTrials - 3));

    trialNames = {};
    for i=1:numTrials
        if ~strcmp(trials(i).name, '.') && ~strcmp(trials(i).name, '..') && (~isempty(strfind(trials(i).name, 'Trial')))
            trialIdx = strsplit(trials(i).name,'Trial');
            trialIdx = str2num(trialIdx{2});
            trialNames{trialIdx} = trials(i).name;
        end
    end


    
    %%% Remove any trials that did not have videos
    newTrialNames = {};
    iter = 1;
    for kk = 1:numel(trialNames)
        if ~isempty(trialNames{kk})
            newTrialNames{iter} = trialNames{kk};
            iter = iter + 1;
        end
    end
    trialNames = newTrialNames;
    subSeq = subSeq(subSeq <= trialIdx);
    subSeq = subSeq(subSeq <= numel(trialNames));
    vids = cell(numel(subSeq), 1);
    viddirs = vids;
    vidnames = vids;
    for kk=subSeq
        disp(['loading trial ', num2str(kk)])
        disp(['loading ', trialNames{kk}])
        
        if isNLW
            viddirs{kk} = fullfile(imageDir, trialNames{kk}, [dataName, '_', trialNames{kk}] );
            
            if ~isempty(zplane)
                viddirs{kk} = fullfile(viddirs{kk}, ['z', num2str(zplane)]);
            end
        else
            viddirs{kk} = fullfile(imageDir, trialNames{kk}, [dataName, '-', trialNames{kk}(end-2:end)] );
        end        
        
        try
%             vidnames{kk} = fullfile(viddirs{kk}, 'reg', 'registered.tif');
            vidnames{kk} = fullfile(viddirs{kk}, 'reg', 'registered_noDenoise.tif');
            vids{kk} = LoadImageStackMultipage(vidnames{kk});
        catch
%             vidnames{kk} = fullfile(viddirs{kk}, 'reg', 'registered_noDenoise.tif');
%             vids{kk} = LoadImageStackMultipage(vidnames{kk});
            error(['No denoised version available for ', vidnames{kk}]);
        end
%         [vids{kk}, vidnames{kk}] = LoadImageStack( viddirs{kk}, [], doLoad );
    end
end
