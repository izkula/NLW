function SaveRegisteredStatsData2p(dataName, imagePath, varargin)
% Save out the max, median, variance, and skewness images across trials
% into the processed data folder. The images are taken as i.e. the maximum
% across trials of the variance within each trial. 

% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('zplanes', [], @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
zplanes = p.Results.zplanes;

%%% Load data
if ~isempty(numTrials)
    subSeq = 1:numTrials
else
    subSeq = []
end

nzplanes = 1;
if ~isempty(zplanes)
    nzplanes = numel(zplanes)
end

for zz = 1:nzplanes
    if ~isempty(zplanes)
        zplane = zplanes(zz);
    else
        zplane = [];
    end
    
    %%% Load all registered videos
    tic
    [vids, vidnames, viddirs] = LoadRegistered2pImageTrials(imagePath, dataName, ...
                                                            'subSeq', subSeq, ...
                                                            'doLoad', true, ...
                                                            'isNLW', isNLW, ...
                                                            'zplane', zplane);
    toc

    imageDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath);


    %%% If you want to save out all trials into a single hdf5 file (doesn't
    %%% actually load much faster than loading directly from tifs, though). 
    % tic
    % h5fname = fullfile(imageDir, 'allRegisteredTrials.h5')
    % H5SaveRegistered2p(h5fname, vids);
    % toc


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
    
    if ~isempty(zplane)
        suffix = ['_z', num2str(zplane), '_'];
    else
        suffix = '';
    end
    
    SaveTiffStack(maxProj, imageDir, ['/maxProj', suffix]);
    SaveTiffStack(maxMean, imageDir, ['/maxMean', suffix]);
    SaveTiffStack(maxVar, imageDir,  ['/maxVar', suffix]);
    SaveTiffStack(maxSkew, imageDir, ['/maxSkew', suffix]);
    save(fullfile(imageDir, ['summaryStatsImages', suffix,'.mat']), 'maxProj', 'maxVar', 'maxSkew', 'maxMean', 'skews', 'vars', 'maxProjs', 'means')

end