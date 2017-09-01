function SaveRegisteredStatsTrials2p(trialPaths, outputPath, varargin)
% Save out the max, median, variance, and skewness images across trials
% into the processed data folder. The images are taken as i.e. the maximum
% across trials of the variance within each trial. 
%%% This function is for generic trial data (i.e. not Bpod). 

% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.parse(varargin{:});
numTrials = p.Results.numTrials;

%%% Load data
if ~isempty(numTrials)
    subSeq = 1:numTrials
else
    subSeq = []
end

%%% Load all registered videos
vids = {};
for i = 1:numel(trialPaths)
    trialPath = trialPaths{i};
    
    if strfind(trialPath, '.tif')
        vids{i} = LoadImageStackMultipage(trialPath, subSeq);
    else
        vids{i} =  LoadImageStack(trialPath, subSeq);
    end        
end


nTrials = numel(vids);
maxProjs = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);
means = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);
vars = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);
skews = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);
for jj = 1:nTrials
    disp(['processing trial ', num2str(jj)])
    vid = vids{jj};
    maxProjs(:,:,jj) = max(vid, [], 3);
    means(:,:,jj) = mean(vid, 3);
    vars(:,:,jj) = var(vid, [], 3);
    skews(:,:,jj) = mean((vid - repmat(means(:,:,jj), 1, 1, size(vid,3))).^3, 3) ...
                    ./( vars(:,:,jj).^(3/2) );
end

maxProj = max(maxProjs, [], 3);
maxMean = max(means, [], 3);
maxVar = max(vars, [], 3);
maxSkew = medfilt2(max(skews, [], 3), [3, 3]);
SaveTiffStack(maxProj, outputPath, '/maxProj');
SaveTiffStack(maxMean, outputPath, '/maxMean');
SaveTiffStack(maxVar, outputPath, '/maxVar');
SaveTiffStack(maxSkew, outputPath, '/maxSkew');
save(fullfile(outputPath, 'summaryStatsImages.mat'), 'maxProj', 'maxVar', 'maxSkew', 'maxMean', 'skews', 'vars', 'maxProjs', 'means')
