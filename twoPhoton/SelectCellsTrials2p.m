function SelectCellsTrials2p(processedDir, varargin)
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
% processedDir: path to folder containing summary statistic images
%%%

global basePath bpodImagePath imagejPath imagejMacroPath

%%% Function options
p = inputParser();
p.addParameter('automateSelection', [], @islogical); 

p.parse(varargin{:});

automateSelection = p.Results.automateSelection;


%%% Load statistical summaries into matlab
S = load(fullfile(processedDir, 'summaryStatsImages.mat'));

if ~automateSelection          

    %%% Do the manual cell selection
    cellPositions = select_cells(S.maxProj, S.maxSkew);
    cellCentroids = cell2mat(cellPositions');
    

    %%% Generate labeled_ROIs
    labeledROIs = zeros(size(S.maxProj));
    labeledROIs = GenerateCellMasksFromCentroids(cellCentroids, 8, 1:size(cellCentroids,1), labeledROIs);
    figure, imagesc(labeledROIs)
    
else
    %%% Do automated segmentation followed by manual adjustment
    thresh = -1;
    [labeledROIs, cellCentroids] = AutoSegmentCells(S.maxProj, S.maxVar, S.maxSkew, thresh);
    [labeledROIs, cellCentroids] = AdjustSegmentation(labeledROIs, cellCentroids, S.maxProj, S.maxSkew);    
end

%%% Plot cell outlines
f= figure, OverlayBoundaries(S.maxProj, labeledROIs);     
% title(strrep(dataName, '_', '-'));
print(fullfile(processedDir, 'cellOutlines.png'), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellOutlines.pdf'));

%%% Plot cell labels
f = figure, imshow(overlayTargets(S.maxProj, cellCentroids), []); colormap jet
% title(strrep(dataName, '_', '-'));
print(fullfile(processedDir, 'cellLabels.png'), '-dpng'); % export_fig(f, fullfile(processedDir, 'cellLabels.pdf'));

%%% Save out cell locations and centroid coordinates
maxProj = S.maxProj;
save(fullfile(processedDir, 'cellPositions.mat'), 'cellCentroids', 'labeledROIs', 'maxProj'); %%% cellCentroids is a vector, each row is a coordinate. labeledROIs is an image. maxProj is an image. 



