clear; clc; close all
init_samz
%
f = dir(fullfile(basePath,'*','*'));

%%
 j = 1

for i = 1:numel(f)
    try
        %get concatenated file names
        if contains(f(i).name,'m876_running')
            filepaths{j} = fullfile(f(i).folder, f(i).name);
            folderpaths{j} = f(i).folder;
            

            j = j+1;
        end
    end
end
%%

 %Define and reset important variables of metadata

datasets = {}; 
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
nidaqSampelRate = 5000;
try
[t,d] = LoadNidaqOutput(fullfile(folderpaths{ii}, 'nidaq.mat'),5);

%extract trial starts
plot(d(2,:))
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 
meta.trialStart = trialTicks(1);
meta.trialEnd = trialTicks(end);

%extract run speed
plot(d(5,:))
runTicks = d(5,:) < -0.1;
tickBin = nidaqSampelRate / 20; %for a bin of 500 ms

runCounts = sum(reshape(runTicks, tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%save out to meta struct
meta.raw.run = d(2,:)
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampelRate
catch
    fprintf(['Cant load nidaq file' name]);
end


if ~exist(fullfile('~', '2presults', 'KX'))
    mkdir(fullfile('~', '2presults', 'KX'))
end

if ~exist(fullfile('~', '2presults', 'KX', [name '_KX']))
    mkdir(fullfile('~', '2presults', 'KX', [name '_KX']))
end
save(fullfile('~', '2presults', 'KX', [name '_KX'],'metadata.mat'),'meta')



    
end
