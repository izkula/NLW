clear; clc; close all
init_samz
%
doRunning = 0;
do8armStimCortex = 0 ;
doGRIN = 0;
do8armFP = 1;

%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doRunning
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

    if~exist(fullfile('~', '2presults', 'Running', aa{end}))
        mkdir(fullfile('~', '2presults', 'Running', aa{end}))
    end
    save(fullfile('~', '2presults', 'Running', aa{end} ,'metadata.mat'),'meta')
end
end
%%%%%%%%%%%%% 8ARM STIM CLAUSTRUM IMAGE CORTEX CODE %%%%%%%%%%%%%%%%
if do8armStimCortex
f = dir(fullfile(basePath,'20180507','*','*'));

j = 1
%get concatenated file names from the session directory
for i = 1:numel(f)
    try
        if contains(f(i).name,'8armstim') && ~contains(f(i).name, 'eye')
            fnames{j} = f(i).name;
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
slice_counter  = 0 %number of slices in the whole image stack
frames = []; %frame numbers for trial starts and ends

load(filepaths{ii})
load(fullfile(folderpaths{ii}, 'AVG_z0_neuron.mat'))

aa = split(folderpaths{ii},'/');
bb = split(aa{end},'_');
name = bb{1};

%get the number of slices in the processed, registed tiff stack
I = imfinfo(fullfile(folderpaths{ii},'AVG_z0_processed.tif'));
nSlices = numel(I);

%load NIDAQ information
try
[t,d] = LoadNidaqOutput(fullfile(folderpaths{ii}, 'nidaq.mat'),5);



                %% CROPPING NIDAQ DATA %%
%Basically because we the imaging screwed up, we are going to end the data
%at the 8th session, yeilding 4 stim and 4 no stim sessions. 
sessionFrames = nSlices*downSampleRate;
lastNidaqTickUsed = sessionFrames / frameSampleRate * nidaqSampleRate;
t = t(1:lastNidaqTickUsed);
d = d(:,1:lastNidaqTickUsed);
frames = info.frame(1:(8*8*2)); %from first to trials tx sessions x on off

                      %% STARTS %%
%extract trial starts
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 
meta.trialStart = trialTicks(1); %trial start
meta.frameStart = frames(1) - round(meta.trialStart / nidaqSampleRate * frameSampleRate) %framestart

                    %% ENDS %% 
meta.trialEnd = trialTicks(end);
meta.frameEnd = meta.frameStart + round(length(d) / nidaqSampleRate * frameSampleRate)


                    %% EXTRACT SPEED %% 
figure; plot(d(5,:))
runTicks = d(5,:) < -0.1;
tickBin = nidaqSampleRate / 40; %for a bin of 100 ms
runCounts = nansum(reshape([runTicks nan(1,84)], tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%1D interpolation
runCountsSmoothAligned = resample(runCountsSmooth, sessionFrames, length(runCountsSmooth));

                %% Stim %%
meta.raw.optoStim = d(1,:);

                %% SAVE OUT TO A META STRUCT %%
meta.raw.run = d(5,:);
meta.run = runCountsSmoothAligned;
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampleRate
meta.frames = frames
catch
    fprintf(['Cant load nidaq file' name]);
end

%     if~exist(fullfile('~', '2presults', '8armStim', aa{end}))
%         mkdir(fullfile('~', '2presults', '8armStim', aa{end}))
%     end
    save(fullfile(folderpaths{ii} ,'metadata.mat'),'meta')
end
end
%%


%%%%%%%%%%%%%%
%%%%%%%%%%%%% GRIN CLAUSTRUM   CODE %%%%%%%%%%%%%%%%
if doGRIN
f = dir(fullfile(basePath,'20180522','*','*'));

downSampleRate = 6;
j = 1
%get concatenated file names from the session directory
for i = 1:numel(f)
    try
        if contains(f(i).name,'8arm2') && ~contains(f(i).name, 'eye')
            fnames{j} = f(i).name;
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
frames = []; %frame numbers for trial starts and ends

load(filepaths{ii})
%load(fullfile(folderpaths{ii}, 'z0_neuron.mat'))

aa = split(folderpaths{ii},'/');
bb = split(aa{end},'_');
name = bb{1};

%get the number of slices in the processed, registed tiff stack
I = imfinfo(fullfile(folderpaths{ii},'AVG_z0_processed.tif'));
nSlices = numel(I);

%load NIDAQ information
try
[t,d] = LoadNidaqOutput(fullfile(folderpaths{ii}, 'nidaq.mat'),5);

try
sessionFrames = nSlices*downSampleRate;
catch
    sessionFrames = nSlices;
end

frames = info.frame; %from first to trials tx sessions x on off

                      %% STARTS %%
%extract trial starts
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 
meta.trialStart = trialTicks(1); %trial start
meta.frameStart = frames(1) - round(meta.trialStart / nidaqSampleRate * frameSampleRate) %framestart

                    %% ENDS %% 
meta.trialEnd = trialTicks(end);
meta.frameEnd = meta.frameStart + round(length(d) / nidaqSampleRate * frameSampleRate)


                    %% EXTRACT SPEED %% 
figure; plot(d(5,:))
runTicks = d(5,:) < -0.1;
tickBin = nidaqSampleRate / 40; %for a bin of 100 ms
runCounts = nansum(reshape(runTicks, tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%1D interpolation
runCountsSmoothAligned = resample(runCountsSmooth, sessionFrames, length(runCountsSmooth));

                %% Stim %%
meta.raw.optoStim = d(1,:);

                %% SAVE OUT TO A META STRUCT %%
meta.raw.run = d(5,:);
meta.run = runCountsSmoothAligned;
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampleRate
meta.frames = frames
catch
    fprintf(['Cant load nidaq file' name]);
    frames = info.frame;
    meta.frameTrial = frames(2) - frames(1);
    meta.frameStart = frames(1) - round(frameSampleRate * 3);
    meta.frameEnd = frames(end) + meta.frameTrial;
    
   
    
    
end

%     if~exist(fullfile('~', '2presults', '8armStim', aa{end}))
%         mkdir(fullfile('~', '2presults', '8armStim', aa{end}))
%     end
    save(fullfile(folderpaths{ii} ,'metadata.mat'),'meta')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 8 ARM FP CLAUSTRUM   CODE %%%%%%%%%%%%%%%%
if do8armFP
f = dir(fullfile(basePath,'20180522','*','*'));

downSampleRate = NaN;
j = 1
%get concatenated file names from the session directory
for i = 1:numel(f)
    try
        if contains(f(i).name,'8arm_fp') && ~contains(f(i).name, 'eye')
            fnames{j} = f(i).name;
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
frames = []; %frame numbers for trial starts and ends

load(filepaths{ii})
%load(fullfile(folderpaths{ii}, 'z0_neuron.mat'))

aa = split(folderpaths{ii},'/');
bb = split(aa{end},'_');
name = bb{1};

%get the number of slices in the processed, registed tiff stack
%I = imfinfo(fullfile(folderpaths{ii},'AVG_z0_processed.tif'));
%nSlices = numel(I);

%load NIDAQ information
try
[t,d] = LoadNidaqOutput(fullfile(folderpaths{ii}, 'nidaq.mat'),5);

try
sessionFrames = 64998;
%catch
%    sessionFrames = nSlices;
end

frames = info.frame; %from first to trials tx sessions x on off

                      %% STARTS %%
%extract trial starts
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 
meta.trialStart = trialTicks(1); %trial start
meta.frameStart = frames(1) - round(meta.trialStart / nidaqSampleRate * frameSampleRate) %framestart

                    %% ENDS %% 
meta.trialEnd = trialTicks(end);
meta.frameEnd = meta.frameStart + round(length(d) / nidaqSampleRate * frameSampleRate)


                    %% EXTRACT SPEED %% 
figure; plot(d(5,:))
runTicks = d(5,:) < -0.04;
tickBin = nidaqSampleRate / 40; %for a bin of 100 ms
runCounts = nansum(reshape(runTicks, tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%1D interpolation
runCountsSmoothAligned = resample(runCountsSmooth, sessionFrames, length(runCountsSmooth));

                %% Stim %%
meta.raw.optoStim = d(1,:);


                %% Photometry Processing %%

meta.raw.fpSignal = d(3,:);
meta.raw.fpcontrol = d(4,:);

meta.fp = d(3,:) - d(4,:);

                

                %% SAVE OUT TO A META STRUCT %%
meta.raw.run = d(5,:);
meta.run = runCountsSmoothAligned;
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampleRate
meta.frames = frames
catch
    fprintf(['Cant load nidaq file' name]);
    %frames = info.frame;
    %meta.frameTrial = frames(2) - frames(1);
    %meta.frameStart = frames(1) - round(frameSampleRate * 3);
    %meta.frameEnd = frames(end) + meta.frameTrial;
    
  
    
end

%     if~exist(fullfile('~', '2presults', '8armStim', aa{end}))
%         mkdir(fullfile('~', '2presults', '8armStim', aa{end}))
%     end
    save(fullfile(folderpaths{ii} ,'metadata.mat'),'meta')
end
end
