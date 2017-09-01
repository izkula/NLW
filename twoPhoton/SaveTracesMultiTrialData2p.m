function SaveTracesMultiTrialData2p(dataName, imagePath, varargin)
% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for one session, containing traces for each trial from the selected cells.
% Traces mat is formatted as: [cell, timepoint, trial]
% 
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath resultsPath

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('zplanes', [], @isnumeric);
p.addParameter('isContinuousTrial',false,@islogical);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
zplanes = p.Results.zplanes;
isContinuousTrial = p.Results.isContinuousTrial;

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


    %%% Load videos
    % [vids, vidnames, viddirs] = Load2pImageTrials(imagePath, dataName, ...
    %                                    'subSeq', subSeq, 'doLoad', true);
    %subSeq = 1:5;
    if ~isContinuousTrial

    [vids, vidnames, viddirs] = LoadRegistered2pImageTrials(imagePath, dataName, ...
                                                            'subSeq', subSeq, ...
                                                            'doLoad', true, ...
                                                            'zplane', zplane);
    else
        viddirs = {fullfile(imagePath, dataName, ['z', num2str(zplane)])};
        error('Not yet implemented')
    end

    metaDatas = {};
    startFrames = zeros(numel(viddirs), 1); %%% Frame at which the signal from bpod arrives (there is only one of these even if imaging across planes)
    useTrial = ones(numel(viddirs), 1);
    nZPlanesRecorded = 1;  %%% Will need to adjust the start frame index accordingly
    if isNLW
        %%% Load imaging metadata file
        for i = 1:numel(viddirs)
            viddir = viddirs{i};
            if ~isempty(zplanes)
                viddir = viddir(1:end-3);
            end
            vidname = strsplit(viddir, '/');
            metaName = fullfile(viddir, [vidname{end}, '__001.mat']);
            if exist(metaName, 'file')
                metaFile = load(metaName);
            else
                warning(['Metadata file does not exist: ', metaName]);
            end
            metaDatas{i} = metaFile.info;
            if metaFile.info.volscan
                nZPlanesRecorded = numel(metaFile.info.otwave_um);
            end
            if isfield(metaFile.info, 'frame')
                startFrames(i) = floor(metaFile.info.frame/nZPlanesRecorded);
                useTrial(i) = 1;
            else
                startFrames(i) = 0;
                useTrial(i) = 0; %%% Applies only to NLW traces
            end
        end
    end
    disp(['nZPlanesRecord: ', num2str(nZPlanesRecorded)]);
    
    if ~isempty(zplane)
        suffix = ['_z', num2str(zplane), '_'];
    else
        suffix = '';
    end


    %%% Load cell coordinates
    % regVidDir = GetProcessedDataDir2p(fullfile(viddirs{1}, 'reg'), basePath, bpodImagePath);
    % load(fullfile(regVidDir, 'cellPositions.mat'), 'cellPositions', 'maxProj')
    processedDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath);
    S = load(fullfile(processedDir, ['cellPositions', suffix, '.mat']));

    saveDir = strrep(imagePath, bpodImagePath, fullfile(resultsPath, '2p/')); %% i.e. ~/Dropbox/oeg_results/2p/gad2m11/Image2POlfGoNoGo/Dec16_2015/Session1
    mkdir(saveDir);

    % nT = 248;
    % nT = size(vids{1}, 3) - 10;
    nTrials = numel(vids)-1;
    nT = min(cellfun(@(a) size(a, 3), vids(1:nTrials)) - startFrames(1:nTrials)); %%% Don't include the last video since it is likely cut off. 
    nCells = max(max(S.labeledROIs)); %size(S.cellCentroids,1);

    traces = zeros(nCells, nT, nTrials);
    dfoftraces = zeros(size(traces));
    fullField = zeros(nT, nTrials);
    dfofFullField = zeros(size(fullField));
    % cellPositions = cell2mat(cellPositions'); % Each row is coordinates of a cell

    maxProjs = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);
    meanProjs = zeros(size(vids{1},1), size(vids{1}, 2), nTrials);

    % compute neuropil regions 
    % construct full activity matrix
    allVids = ConcatenateVids(vids(1:nTrials), nT, startFrames);
    doRemoveOutliers = false;
    [neuropilTraces,neuropilMasks] = FindNeuropilTraces(S.labeledROIs, allVids, vids, ...
                                                        nT, nTrials, startFrames, doRemoveOutliers);
    alpha = 0.7;

    for jj = 1:nTrials
        vid = vids{jj};
        if size(vids{jj}, 3) < nT
            break
        end
        maxProjs(:,:,jj) = max(vid, [], 3);
        meanProjs(:,:,jj) = mean(vid, 3);

        tt = ExtractTraces2p(S.labeledROIs, vid);
        traces(:,:,jj) = tt(:, 1+startFrames(jj):nT+startFrames(jj));

        %%% Do deltaF of the traces
        M = repmat(median(traces(:,:,jj), 2), 1, size(traces,2));
        dfoftraces(:,:,jj) = (traces(:,:, jj) -M)./M;

        fullField(:, jj) = squeeze(sum(squeeze(sum(vid(:,:,1+startFrames(jj):nT+startFrames(jj)), 1)), 1));
        Mf = median(fullField(:,jj));
        dfofFullField(:,jj) = (fullField(:, jj) - Mf)./Mf;
        %figure, plot(traces(:,:,jj)'); 
    end

    maxProj = max(maxProjs, [], 3);
    meanProj = mean(meanProjs, 3);
    S.maxProj = maxProj;
    S.meanProj = meanProj;
    %%% Plot overlay of cell positions with max project of that trial
    f = figure, imshow(overlayTargets(S.maxProj, S.cellCentroids), []); colormap jet
    export_fig(f, fullfile(saveDir, ['cellPositionsMax', suffix, '.pdf']));
    f = figure, imshow(overlayTargets(S.meanProj, S.cellCentroids), []); colormap jet
    export_fig(f, fullfile(saveDir, ['cellPositionsMean', suffix, '.pdf']));

    %%% Save out traces matrix
    cellCentroids = S.cellCentroids;
    labeledROIs = S.labeledROIs;
    Ncells = size(cellCentroids, 1);
    mkdir(saveDir);
    save(fullfile(saveDir, ['traces', suffix, '.mat']), 'neuropilTraces', 'neuropilMasks', ...
                                           'traces', 'dfoftraces', 'dfofFullField', ...
                                           'cellCentroids', 'labeledROIs', 'maxProj', ...
                                           'meanProj', 'fullField', 'Ncells', ...
                                           'metaDatas', 'useTrial');
end                                   
                                   


                                   
                                   
function allVids = ConcatenateVids(v, nT, startFrames)
%%% v - cell array of vids to concatenate
%%% nT - number of time points to uniformly trim videos
%%% startFrames - frame corresponding to start of bpod acquisition

temp = v{1};
% nT = size(temp,3)-2;
allVids = zeros(size(temp,1), size(temp,2), nT*numel(v));
tic
k = 1;
for i=1:numel(v)
    i
    if size(v{i}, 3) < nT
        break
    end
    allVids(:,:,(i-1)*nT+1:i*nT) = v{i}(:,:,1+startFrames(i):nT+startFrames(i));
end
toc


function [neuropilTraces, neuropilMasks] = FindNeuropilTraces(cellMasks, ...
                                        allVids, vids, nT, nTrials, startFrames, ...
                                        doRemoveOutliers)
disp('calculating neuropil traces');
% takes concated vids and img mask
Ncells = max(cellMasks(:));

masks = cellMasks;
[Nx,Ny] = size(cellMasks);
traces = ExtractTraces2p(cellMasks, allVids);
neuropilTraces = cell(numel(vids),1);
neuropilMasks = zeros(size(cellMasks));
cellInd = cell(Ncells,1);
if Nx == 512
    strelsize = 6;
else
    strelsize = 4;
end
for i=1:Ncells
    ind = find(cellMasks == i);
    currMask = zeros(Nx,Ny);
    %[x1,x2] = ind2sub(size(cellMasks), inds(i));
    % construct donut mask
    currMask(ind) = 1; 
    se = strel('disk',strelsize);
    dilated = imdilate(currMask,se);
    dilated(ind) = 0;
    dilated(find(cellMasks > 0)) = 0; %%% Remove any neighboring cells
    donutInd = find(dilated == 1);
    currTrace = traces(i,:);
    
    % remove contaminating pixels (i.e. any pixels in the mask that have
    % traces that fluctuate too extensively (i.e. that likely overlap with 
    % neuron soma).
    if doRemoveOutliers
        for j=1:numel(donutInd)
            [x,y] = ind2sub(size(cellMasks), donutInd(j));
            indTrace = allVids(x,y,:);
            deltaTrace = indTrace(:) - currTrace(:);
            traceSD = std(deltaTrace(:));
            if any(deltaTrace > 2*traceSD)
                dilated(x,y) = 0;
            end                
        end
    end
    ind = find(dilated==1);
    neuropilMasks(dilated==1) = i;
    cellInd{i} = ind;
end
for i=1:nTrials  
    currVid = vids{i};
    currVid = reshape(currVid, [size(currVid,1) * size(currVid,2), size(currVid,3)]);
%     traces = zeros(Ncells, size(vids{i},3));
    traces = zeros(Ncells, nT);
    for j=1:Ncells        
        ind = cellInd{j};
        traces(j,:) = squeeze(mean(currVid(ind,1+startFrames(i):nT+startFrames(i)),1));
    end
    neuropilTraces{i} = traces;
end
