function MakeAverageVid(viddir, vidfile, outdir, whichTrials)
%%% Make average video
%%% viddir - directory that also contains the trialStarts mat file
%%% vidfile - the rest of the path from viddir down to the vid file,
%%% including .tif in the name
%%% outdir - where to save the movie
%%% whichTrials - if only want to average across a subset

if ~exist('whichTrials', 'var') || isempty(whichTrials)
   doAllTrials = true; 
else
   doAllTrials = false;
end

% viddir = '/media/Data/processedData/BpodImageData/ckrspm2/ImageNLWTrainStimGoNoGo/Jun01_2017/Session1/ckrspm2_ImageNLWTrainStimGoNoGo_Jun01_2017_Session1/'
% viddir = '/media/Data/processedData/BpodImageData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session1/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session1/'
% viddir = '/media/Data/processedData/BpodImageData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session2/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session2/'

vid = LoadImageStackMultipage(fullfile(viddir, 'z0', 'reg', 'registered_noDenoise.tif'));

TS = load(fullfile(viddir, 'z0_trialStarts.mat'));

%%

startFrames = TS.trialStartInds;
endFrames = TS.trialEndInds;
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

s = round(startFrames/3) + 1;
e = round(endFrames/3) + 1;

% s = round(TS.trialStartInds/3)+1;
% e = round(TS.trialEndInds/3)+1;

v = zeros(size(vid, 1), size(vid, 2), max(e - s)+1);

if doAllTrials
    whichTrials = [1:numel(e)-2];
else
   whichTrials = whichTrials(1:end-1); 
end
   

trial_vids = zeros(size(v, 1), size(v, 2), size(v, 3), numel(whichTrials));
% for i = whichTrials
for kk = 1:numel(whichTrials)
    i = whichTrials(kk);
    v(:,:,1:(e(i)-s(i)+1)) = v(:,:,1:(e(i)-s(i)+1)) + vid(:,:,s(i):e(i));
    trial_vids(:,:,1:(e(i)-s(i)+1), kk) = vid(:,:,s(i):e(i));
end

v = v/numel(whichTrials);

%%
% fnames = SaveTiffStack( v, fullfile(viddir, 'trial_average'))
fnames = SaveTiffStack( v, fullfile(outdir, 'trial_average'))