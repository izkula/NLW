function SaveTracesTrials2p(trialPaths, processedDir, saveDir, varargin)
% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for one session, containing traces for each trial from the selected cells.
% Traces mat is formatted as: [cell, timepoint, trial]
% 
% Parameters:
% trialPaths: cell array with paths to each video 
% processedDir: path to folder containing cell locations file
% saveDir: output path for saved traces (i.e. dropbox)
%%%

global basePath bpodImagePath resultsPath

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
% p.addParameter('cellRadius', 3, @isnumeric); % In pixels. Set to [] if using automatic cell selection. 
p.parse(varargin{:});
numTrials = p.Results.numTrials;
% cellRadius = p.Results.cellRadius;

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

                                                    
                                                    
    
                                                    
%%% Load cell coordinates
S = load(fullfile(processedDir, 'cellPositions.mat'));

mkdir(saveDir);

nT = size(vids{1}, 3);
nTrials = numel(vids);
nCells = max(max(S.labeledROIs)); %size(S.cellCentroids,1);

traces = zeros(nCells, nT, nTrials);
dfoftraces = zeros(size(traces));
fullField = zeros(nT, nTrials);
dfofFullField = zeros(size(fullField));
% cellPositions = cell2mat(cellPositions'); % Each row is coordinates of a cell

maxProjs = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);

% compute neuropil regions 
% construct full activity matrix
allVids = ConcatenateVids(vids);
[neuropilTraces,neuropilMasks] = FindNeuropilTraces(S.labeledROIs, allVids, vids);
alpha = 0.7;

for jj = 1:nTrials
    vid = vids{jj};
    maxProjs(:,:,jj) = max(vid, [], 3);
    %%% Plot overlay of cell positions with max project of that trial
%     close all
%     figure, imagesc(max(vid, [], 3)); hold on; plot(cellPositions(:,1), cellPositions(:,2), 'r+');
    
    traces(:,:,jj) = ExtractTraces2p(S.labeledROIs, vid);
    
    %%% Do deltaF of the traces
    M = repmat(median(traces(:,:,jj), 2), 1, size(traces,2));
    dfoftraces(:,:,jj) = (traces(:,:, jj) -M)./M;
    
    fullField(:, jj) = squeeze(sum(squeeze(sum(vid, 1)), 1));
    Mf = median(fullField(:,jj));
    dfofFullField(:,jj) = (fullField(:, jj) - Mf)./Mf;
    %figure, plot(traces(:,:,jj)'); 
end

maxProj = max(maxProjs, [], 3);
f = figure, imshow(overlayTargets(S.maxProj, S.cellCentroids), []); colormap jet
print(fullfile(saveDir, 'cellPositions.png'), '-dpng')
export_fig(f, fullfile(saveDir, 'cellPositions.pdf'));

%%% Save out traces matrix
cellCentroids = S.cellCentroids;
labeledROIs = S.labeledROIs;
Ncells = size(cellCentroids, 1);
save(fullfile(saveDir, 'traces.mat'), 'neuropilTraces', 'neuropilMasks', 'traces', 'dfoftraces', 'dfofFullField', 'cellCentroids', 'labeledROIs', 'maxProj', 'fullField', 'Ncells');

function allVids = ConcatenateVids(v)
temp = v{1};
nT = size(temp,3);
allVids = zeros(size(temp,1), size(temp,2), nT*numel(v));
k = 1;
for i=1:numel(v)
    for j=1:nT
        allVids(:,:,k) = v{i}(:,:,j);
        k = k + 1;
    end
end

function [neuropilTraces, neuropilMasks] = FindNeuropilTraces(cellMasks, allVids, vids)
disp('calculating neuropil traces');
% takes concated vids and img mask
Ncells = max(cellMasks(:));

masks = cellMasks;
[Nx,Ny] = size(cellMasks);
traces = ExtractTraces2p(cellMasks, allVids);
neuropilTraces = cell(numel(vids),1);
neuropilMasks = zeros(size(cellMasks));
cellInd = cell(Ncells,1);
for i=1:Ncells
    ind = find(cellMasks == i);
    currMask = zeros(Nx,Ny);
    %[x1,x2] = ind2sub(size(cellMasks), inds(i));
    % construct donut mask
    currMask(ind) = 1; 
    se = strel('disk',4);
    dilated = imdilate(currMask,se);
    dilated(ind) = 0;
    donutInd = find(dilated == 1);
    currTrace = traces(i,:);
    % remove contaminating traces
    for j=1:numel(donutInd)
        [x,y] = ind2sub(size(cellMasks), donutInd(j));
        indTrace = allVids(x,y,:);
        deltaTrace = indTrace(:) - currTrace(:);
        traceSD = std(deltaTrace(:));
        if any(deltaTrace > 2*traceSD)
            dilated(x,y) = 0;
        end                
    end
    ind = find(dilated==1);
    neuropilMasks(dilated==1) = i;
    cellInd{i} = ind;
end
for i=1:numel(vids)   
    currVid = vids{i};
    currVid = reshape(currVid, [size(currVid,1) * size(currVid,2), size(currVid,3)]);
    traces = zeros(Ncells, size(vids{i},3));
    for j=1:Ncells        
        ind = cellInd{j};
        traces(j,:) = squeeze(mean(currVid(ind,:),1));
    end
    neuropilTraces{i} = traces;
end
