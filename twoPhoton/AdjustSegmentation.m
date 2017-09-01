function [labeledROIs, newCellCentroids] = AdjustSegmentation(labeledROIs, ...
                                           cellCentroids, f_max, f_skew, varargin)
%%% Provided the output of AutoSegmentCells, manually adjust cell
%%% selection. Return modified labeledROIs and cellCentroids (same format
%%% as from AutoSegmentCells).
       
p = inputParser();
p.addParameter('cellRadius', [], @isnumeric);

p.parse(varargin{:});
cellRadius = p.Results.cellRadius;

% close all
newCellCentroids = cell2mat(select_cells(f_max, f_skew, num2cell(cellCentroids,2), labeledROIs)');
[new_points, removed_points, new_inds, removed_inds] = GetPointChanges(num2cell(newCellCentroids,2), num2cell(cellCentroids,2));
%%% remove those points from the labeled image
for i = 1:numel(removed_inds)
    labeledROIs(labeledROIs == removed_inds(i)) = 0; 
end
%%% add new points to the labeled image
if numel(new_inds) > 0
    labeledROIs = GenerateCellMasksFromCentroids(new_points, cellRadius, new_inds, labeledROIs);
end
