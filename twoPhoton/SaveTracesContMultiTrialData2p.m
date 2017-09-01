function SaveTracesContMultiTrialData2p(dataName, imagePath, varargin)
% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for one session, containing traces for each trial from the selected cells.
% This function is for trials that were recorded continuously (as opposed
% to splitting up each trial into a separate vid). 
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
p.addParameter('isContinuousTrial',true,@islogical);
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
zplanes = p.Results.zplanes;
isContinuousTrial = p.Results.isContinuousTrial;
tDownsampleFactor = p.Results.tDownsampleFactor;

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

tDownsampleFactor = tDownsampleFactor/nzplanes;

for zz = 1:nzplanes
    if ~isempty(zplanes)
        zplane = zplanes(zz);
    else
        zplane = [];
    end

    if ~isempty(zplane)
        zstr = ['z', num2str(zplane)];
    else
        zstr = '';
    end


    processedDataDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath)

    viddir = fullfile(processedDataDir, dataName, zstr);
    vid = LoadImageStackMultipage(fullfile(viddir, 'reg', 'registered_noDenoise.tif'));

    
% % % % Account for the fact that you temporally downsampled. 
% % % % Extract frame time for each trial start.         
        
    %%% Extract trial start times (accounting for the temporal downsampling)
    if exist(fullfile(processedDataDir, dataName, [zstr, '_trialStarts.mat']), 'file')
        TS = load(fullfile(processedDataDir, dataName, [zstr, '_trialStarts.mat']))
        startFrames = TS.trialStartInds; %%% Does not account for temporaldownsampling
        endFrames = TS.trialEndInds;
        useTrial = TS.didNLWRecordStartFrame;
        metaDatas = {};
        metaDatas{1}.volscan = true; %%% Temporary hack for a workaround for future functions. 
        metaDatas{1}.otwave_um = zeros(1, nzplanes);
    else
       error('Not yet implemented') 
%         metaDatas = {};
%         startFrames = zeros(numel(viddirs), 1); %%% Frame at which the signal from bpod arrives (there is only one of these even if imaging across planes)
%         useTrial = ones(numel(viddirs), 1);
%         nZPlanesRecorded = 1;  %%% Will need to adjust the start frame index accordingly
%         if isNLW
%             %%% Load imaging metadata file
%             for i = 1:numel(viddirs)
%                 viddir = viddirs{i};
%                 if ~isempty(zplanes)
%                     viddir = viddir(1:end-3);
%                 end
%                 vidname = strsplit(viddir, '/');
%                 metaName = fullfile(viddir, [vidname{end}, '__001.mat']);
%                 if exist(metaName, 'file')
%                     metaFile = load(metaName);
%                 else
%                     warning(['Metadata file does not exist: ', metaName]);
%                 end
%                 metaDatas{i} = metaFile.info;
%                 if metaFile.info.volscan
%                     nZPlanesRecorded = numel(metaFile.info.otwave_um);
%                 end
%                 if isfield(metaFile.info, 'frame')
%                     startFrames(i) = floor(metaFile.info.frame/nZPlanesRecorded);
%                     useTrial(i) = 1;
%                 else
%                     startFrames(i) = 0;
%                     useTrial(i) = 0; %%% Applies only to NLW traces
%                 end
%             end
%         end
%         disp(['nZPlanesRecord: ', num2str(nZPlanesRecorded)]);
    end
        
    %%% Load cell coordinates
    if ~isempty(zplane)
        suffix = ['_z', num2str(zplane), '_'];
    else
        suffix = '';
    end
    S = load(fullfile(processedDataDir, dataName, zstr, 'reg', ['cellPositions', suffix, '.mat']));

%     saveDir = strrep(imagePath, bpodImagePath, fullfile(resultsPath, '2p/')); %% i.e. ~/Dropbox/oeg_results/2p/gad2m11/Image2POlfGoNoGo/Dec16_2015/Session1
    saveDir = strrep(imagePath, bpodImagePath, fullfile(resultsPath, '2p')); %% i.e. ~/Dropbox/oeg_results/2p/gad2m11/Image2POlfGoNoGo/Dec16_2015/Session1
    mkdir(saveDir);
    
    
    %%% Extract traces
    nCells = max(max(S.labeledROIs)); %size(S.cellCentroids,1);
    nT = size(vid, 3);

    traces = zeros(nCells, nT);
    dfoftraces = zeros(size(traces));
    fullField = zeros(nT);
    dfofFullField = zeros(size(fullField));

    maxProjs = zeros(size(vid,1), size(vid, 2));
    meanProjs = zeros(size(vid,1), size(vid, 2));
    
    maxProj = max(vid, [], 3);
    meanProj = mean(vid, 3);

    traces = ExtractTraces2p(S.labeledROIs, vid);

    %%% Extract neuropil
    doRemoveOutliers = false;
    [neuropilTraces,neuropilMasks] = FindNeuropilTraces(S.labeledROIs, traces, vid, ...
                                                        nT, doRemoveOutliers);
    
    %%% Do deltaF of the traces
    M = median(traces, 2);
    dfoftraces = (traces -M)./M;

    M = median(neuropilTraces, 2);
    dfofNeuropilTraces = (neuropilTraces -M)./M;
    
    fullField = squeeze(sum(squeeze(sum(vid, 1)), 1));
    Mf = median(fullField);
    dfofFullField = (fullField - Mf)./Mf;    
    
    
    %%% Do some reorganization for compatibility with later functions
    appendedData = {}; %%% All traces appended (as originally recorded the data). 
    appendedData.appendedTraces = traces;
    appendedData.appendedNeuropil = neuropilTraces;
    appendedData.appendedFullFild = fullField;
    
    
    
%     %%%% Check for overflow error
    if max(abs(diff(startFrames))) > 20000
        [~, ind] = max(abs(startFrames));
        startFrames_new = startFrames;
        startFrames_new(ind+1:end) = startFrames_new(ind+1:end) + double(intmax('uint16'));
        
        [~, ind] = max(abs(endFrames));
        endFrames_new = endFrames;
        endFrames_new(ind+1:end) = endFrames_new(ind+1:end) + double(intmax('uint16'));
        
        startFrames = startFrames_new;
        endFrames = endFrames_new;
    end

%     if max(abs(diff(sf)))> 6000 %%% chose a random number that accounts for some potentially long trials
%         [~, ind] = max(abs(sf))
%         sf_new = sf
%         ef_new = ef;
%         dt_sf = round(mean(diff(sf(ind-30:ind))))
%         sf_new(ind+1:end) = sf_new(ind+1:end) + sf(ind) - sf(ind+1) + dt_sf;
%         ef_new(ind:end) = ef_new(ind:end) + ef(ind-1) - ef(ind) + dt_sf;
%         
%         sf = sf_new;
%         ef = ef_new;
%     end
    
    sf = ceil(startFrames/tDownsampleFactor/nzplanes);
    ef = ceil(endFrames/tDownsampleFactor/nzplanes);
    traces = ReshapeToTrials(traces, sf, ef);
    dfoftraces = ReshapeToTrials(dfoftraces, sf, ef);
    neuropilTraces = ReshapeToTrials(neuropilTraces, sf, ef); %%% For compatibility with later functions
    dfofNeuropilTraces = ReshapeToTrials(dfofNeuropilTraces, sf, ef);
    
    fullField = squeeze(ReshapeToTrials(fullField, sf, ef));
    dfofFullField = squeeze(ReshapeToTrials(dfofFullField, sf, ef));

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
                                           'traces', 'dfoftraces', 'dfofNeuropilTraces', 'dfofFullField', ...
                                           'cellCentroids', 'labeledROIs', 'maxProj', ...
                                           'meanProj', 'fullField', 'Ncells', ...
                                           'metaDatas', 'useTrial', ...
                                           'startFrames', 'endFrames', ...
                                           'appendedData', 'tDownsampleFactor');
end                                   
                                   


function [neuropilTraces, neuropilMasks] = FindNeuropilTraces(cellMasks, ...
                                        traces, vid, nT, ...
                                        doRemoveOutliers)
disp('calculating neuropil traces');
% takes concated vids and img mask
Ncells = max(cellMasks(:));

masks = cellMasks;
[Nx,Ny] = size(cellMasks);

neuropilMasks = zeros(size(cellMasks));
cellInd = cell(Ncells,1);
if Nx == 512
    strelsize = 6;
else
    strelsize = 6;
end
for i=1:Ncells
    ind = find(cellMasks == i);
    currMask = zeros(Nx,Ny);

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
            indTrace = vid(x,y,:);
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

currVid = vid;
currVid = reshape(currVid, [size(currVid,1) * size(currVid,2), size(currVid,3)]);
traces = zeros(Ncells, nT);
for j=1:Ncells        
    ind = cellInd{j};
    traces(j,:) = squeeze(mean(currVid(ind,:),1));
end
neuropilTraces = traces;
