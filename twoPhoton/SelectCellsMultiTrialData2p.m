function SelectCellsMultiTrialData2p(dataName, imagePath, varargin)
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

p.parse(varargin{:});

automateSelection = p.Results.automateSelection;
numTrials = p.Results.numTrials;
cellRadius = p.Results.cellRadius;
zplanes = p.Results.zplanes;
doUseVideo = p.Results.doUseVideo;

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

    if doUseVideo
        [~, ~, viddirs] = Load2pImageTrials(imagePath, dataName, ...
                                           'subSeq', subSeq, ...
                                           'doLoad', false, ...
                                           'isNLW', isNLW, ...
                                           'zplane', zplane);
    end
                                   
    %%% Load statistical summaries into matlab
    processedDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath);
    
    if ~isempty(zplane)
        suffix = ['_z', num2str(zplane), '_'];
    else
        suffix = '';
    end
    
    
    S = load(fullfile(processedDir, ['summaryStatsImages', suffix, '.mat']));

    if isempty(cellRadius)
        if size(S.maxProj, 1) == 512
            cellRadius = 6; %when using the 16x on Bruker
        else
            cellRadius = 4;        
        end
    end

    if ~automateSelection          
        if doUseVideo
             regVidDir = GetProcessedDataDir2p(fullfile(viddirs{1}, 'reg'), basePath, bpodImagePath);
            regVidFile = fullfile(regVidDir, 'registered.tif')
            %%% Load video
            vid = LoadImageStackMultipage(regVidFile);

            %%% Select cells
            maxProj = max(vid, [], 3);
            cellPositions = select_cells(S.maxProj, vid);
        end

        %%% Do the manual cell selection
        cellPositions = select_cells(S.maxMean, S.maxSkew);
        cellCentroids = cell2mat(cellPositions');


        %%% Generate labeled_ROIs
        labeledROIs = zeros(size(S.maxProj));
        labeledROIs = GenerateCellMasksFromCentroids(cellCentroids, cellRadius, 1:size(cellCentroids,1), labeledROIs);
        figure, imagesc(labeledROIs)

    else
        %%% Do automated segmentation followed by manual adjustment
    %     [labeledROIs, cellCentroids] = AutoSegmentCells(S.maxProj, S.maxVar, S.maxSkew);
        [labeledROIs, cellCentroids] = AutoSegmentCells(S.maxMean, S.maxVar, S.maxSkew);
        [labeledROIs, cellCentroids] = AdjustSegmentation(labeledROIs, cellCentroids, S.maxMean, S.maxSkew, 'cellRadius', cellRadius);    
    end
    
    while 1
        figure('Position', [1000, 1000, 1300, 1300]), OverlayBoundaries(S.maxMean, labeledROIs);
        disp(['Num cells: ', num2str(numel(unique(labeledROIs) - 1))])
        str = input('Are you happy with segmentation? (type n for no, y for yes) ','s')
        if strcmp(str, 'n')
           cellRadInput = input(['What cell radius should be used? (current is ', num2str(cellRadius), ') '],'s');
           [labeledROIs, cellCentroids] = AdjustSegmentation(labeledROIs, cellCentroids, S.maxMean, S.maxSkew, 'cellRadius', str2double(cellRadInput));    
        else
            break;
        end
    end

    f= figure, imagesc(labeledROIs);     
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellMasks', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellOutlines.pdf'));

    %%% Plot cell outlines
    f= figure, OverlayBoundaries(S.maxMean, labeledROIs);     
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellOutlines', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellOutlines.pdf'));

    %%% Plot cell labels
    f = figure, imshow(overlayTargets(S.maxMean, cellCentroids), []); colormap jet
    title(strrep(dataName, '_', '-'));
    print(fullfile(processedDir, ['cellLabels', suffix, '.png']), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellLabels.pdf'));

    %%% Save out cell locations and centroid coordinates
    maxProj = S.maxProj;
    save(fullfile(processedDir, ['cellPositions', suffix, '.mat']), 'cellCentroids', 'labeledROIs', 'maxProj'); %%% cellCentroids is a vector, each row is a coordinate. labeledROIs is an image. maxProj is an image. 

end

