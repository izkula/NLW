clear; clc; close all
init_samz
%
doRunning = 1;
do8arm = 0 ;
doCortex = 0;



f = dir(fullfile(basePath,'*','*'));

%%
 j = 1

for i = 1:numel(f)
    try
        %get concatenated file names
        if contains(f(i).name,'running')
            filepaths{j} = fullfile(f(i).folder, f(i).name);
            folderpaths{j} = f(i).folder;
            
            j = j+1;
        end
    end
end
%%

 %Define and reset important variables of metadata

datasets = {}; frameSampleRate = 30.98; nidaqSampleRate = 5000;


for ii = 1:numel(filepaths)
slice_counter = 0 %number of slices in the whole image stack
frames = []; %frame numbers for trial starts and ends

load(filepaths{ii})

aa = split(folderpaths{ii},'/');
bb = split(aa{end},'_');
name = bb{1};

%get the number of slices in the processed, registed tiff stack
I = imfinfo(fullfile(folderpaths{ii},'z0_processed.tif'));
nSlices = numel(I);
try
frames = [frames; info.frame + slice_counter];
plot(frames)
meta.frames = frames;

catch
    fprintf(['no frames recorded for ' name])
end

%add this to slice counter
slice_counter = slice_counter + nSlices


%load NIDAQ information
try
[t,d] = LoadNidaqOutput(fullfile(folderpaths{ii}, 'nidaq.mat'),5);

%extract trial starts
plot(d(2,:))
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 
meta.trialStart = trialTicks(1);
meta.trialEnd = trialTicks(end);


%extract the imaging frame the corresponds to the beginning of the
%experiment
meta.frameStart = meta.frames(1) - round(meta.trialStart / nidaqSampleRate * frameSampleRate)
meta.frameEnd = meta.frames(1) + round(length(d) / nidaqSampleRate * frameSampleRate)


%extract run speed
figure; plot(d(5,:))
runTicks = d(5,:) < -0.1;
tickBin = nidaqSampleRate / 20; %for a bin of 50 ms
runCounts = sum(reshape(runTicks, tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%1D interpolation
sessionFrames = (meta.frameEnd - meta.frameStart);
runCountsSmoothAligned = resample(runCountsSmooth, sessionFrames, length(runCountsSmooth));

%save out to meta struct
meta.raw.run = d(2,:)
meta.run = runCountsSmoothAligned;
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampleRate
catch
    fprintf(['Cant load nidaq file' name]);
end

if doRunning
    if~exist(fullfile('~', '2presults', 'Running', aa{end}))
        mkdir(fullfile('~', '2presults', 'Running', aa{end}))
    end
    save(fullfile('~', '2presults', 'Running', aa{end} ,'metadata.mat'),'meta')
end

if do8arm
end

if doCortex 
end




    
end
