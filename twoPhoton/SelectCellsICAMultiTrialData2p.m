function SelectCellsICAMultiTrialData2p(dataName, imagePath, varargin)
% Provides a gui interface for selecting cells from a single specified
% data session (which consists of multiple trials).
% The assumption is that the videos from trials within the session have
% been registered, and thus cell selection only need occur using the first
% file. 
%
% Saves a mat with cell locations to the corresponding  processedData folder 
% for that session. 
%
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e. path to the raw dataset
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath imagejPath imagejMacroPath

%%% Function options
p = inputParser();
p.addParameter('automateSelection', [], @islogical); 
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('cellRadius', [], @isnumeric);
p.addParameter('zplanes', [], @isnumeric);
p.addParameter('doUseVideo', false, @islogical); %%% Use video from a single trial for help with manual selection
p.addParameter('doAppendTiffs', false, @islogical); %%% Use video from a single trial for help with manual selection
p.addParameter('isNLW', true, @islogical); %%% Use video from a single trial for help with manual selection

p.parse(varargin{:});

automateSelection = p.Results.automateSelection;
numTrials = p.Results.numTrials;
cellRadius = p.Results.cellRadius;
zplanes = p.Results.zplanes;
doUseVideo = p.Results.doUseVideo;
doAppendTiffs = p.Results.doAppendTiffs;
isNLW = p.Results.isNLW;

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
    
    if ~isempty(zplane)
        suffix = ['_z', num2str(zplane), '_'];
    else
        suffix = '';
    end
    
    
    %%% Load ICA based segmentation
    currdir = {fullfile(imagePath, dataName, ['z', num2str(zplane)])};
    processedDir = GetProcessedDataDir2p(fullfile(currdir, 'reg'), basePath, bpodImagePath)
    processedDir = processedDir{1};
    segdir = fullfile(processedDir, 'timecourses');
    segfname = dir(fullfile(segdir, '*_cellmasks.mat'));
    seg = load(fullfile(segdir, segfname.name));

    movm = seg.movm;
    
    labeledROIs = zeros(size(movm));
    for a = 1:size(seg.cellmask, 1)
        labeledROIs(seg.cellmask(a, :,:) > 0) = a;
    end
    cellCentroids = seg.segcentroid;
    figure, imagesc(labeledROIs + 100*(labeledROIs > 0))
    
    %%% Now adjust the segmentation
    %     [labeledROIs, cellCentroids] = AutoSegmentCells(S.maxProj, S.maxVar, S.maxSkew);
    [labeledROIs, cellCentroids] = AdjustSegmentation(labeledROIs, cellCentroids, movm, movm, 'cellRadius', cellRadius);    
    
    while 1
        figure('Position', [1000, 1000, 1300, 1300]), OverlayBoundaries(movm, labeledROIs);
        disp(['Num cells: ', num2str(numel(unique(labeledROIs) - 1))])
        str = input('Are you happy with segmentation? (type n for no, y for yes) ','s')
        if strcmp(str, 'n')
           cellRadInput = input(['What cell radius should be used? (current is ', num2str(cellRadius), ') '],'s');
           [labeledROIs, cellCentroids] = AdjustSegmentation(labeledROIs, cellCentroids, movm, movm, 'cellRadius', str2double(cellRadInput));    
        else
            break;
        end
    end
    

    f= figure, imagesc(labeledROIs);     
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellMasks', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellOutlines.pdf'));

    %%% Plot cell outlines
    f= figure, OverlayBoundaries(movm, labeledROIs);     
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellOutlines', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellOutlines.pdf'));

    %%% Plot cell labels
    f = figure, imshow(overlayTargets(movm, cellCentroids), []); colormap jet
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellLabels', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellLabels.pdf'));

    %%% Save out cell locations and centroid coordinates
%     maxProj = S.maxProj;
    maxProj = movm;
    save(fullfile(processedDir, ['cellPositions', suffix, '.mat']), 'cellCentroids', 'labeledROIs', 'maxProj'); %%% cellCentroids is a vector, each row is a coordinate. labeledROIs is an image. maxProj is an image. 
end

