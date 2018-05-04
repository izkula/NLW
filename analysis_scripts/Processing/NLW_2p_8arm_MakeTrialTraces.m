%%This analysis script makes the trial averaged data for every claustrum
%%neuron, saves out a trialtrace.mat file%%%
clear all; close all; clc



init_samz

f = dir(fullfile('~/2presults/', '8arm'))

trialOrder = {1 2 1 4 2 1 2 3};

for i = 9:numel(f)
    dirName = fullfile(f(i).folder,f(i).name);
    
    %load neuron trace data
    load(fullfile(dirName, 'neuron.mat'))

    %load meta data (frame times)
    load(fullfile(dirName, 'metadata.mat'))

    %Hard Code some things for the 8arm experiment
    subject = dirName(end-8:end-5);
    nNeurons = size(C,1);
    nFramesBeforeTrial = nFramesBaseline;
    nFramesAfterTrial =nFramesPostTrial;
    nFramesPerTrial = nFramesAdvance + nFramesStationary + nFramesRetreat + nFramesBeforeTrial + nFramesAfterTrial;
    nArms = 8;
    
    %one mouse had fewer trials.
    if contains(dirName(end-8:end-5),'m874')
        nTrials = 8
    else
        nTrials = 10
    end
    
        
    
    %extract traces 
    %note, from frame start to frame end, 
        frames = meta.frames;
        
        %sometimes the first frame is 0.
        if frames(1) == 0
            frames = frames(2:end)
        end
        
        %split traces into trials
        traces_preliminary = [];
        
        for j = 1:size(frames,1)/2
            frame_start = frames(2*j -1) - nFramesBeforeTrial;
            frame_end = frames(2*j) + nFramesAfterTrial;
            
            try
           if frame_end-frame_start == nFramesPerTrial
                traces_preliminary(:,:,j) = C(:,frame_start:frame_end); 
    
           elseif frame_end-frame_start == nFramesPerTrial - 1
                traces_preliminary(:,:,j) = C(:,frame_start:frame_end+1); %if off by 1 frame
        
           elseif frame_end - frame_start == nFramesPerTrial + 1
                traces_preliminary(:,:,j) = C(:,frame_start:frame_end-1); %if off by 1 frame

           else
                fprintf(num2str(j))
                printf([' A potentional trial is skipped'])           
           end
            catch 
                fprintf('Error')
            end
            
        end              
    
    
    %organize traces by type of stimulus
    traces = cell(1,nArms);
    trialidx = reshape(repmat(1:nTrials,nArms,1),1,[]);
    idx = repmat(1:nArms,1,nTrials);
    
    for ii = 1:size(traces_preliminary,3)
        traces{idx(ii)}(:,:,trialidx(ii))= traces_preliminary(:,:,ii);
    end
    
    save(fullfile(f(i).folder, f(i).name, 'trialtraces.mat'), ...
        'traces', 'traces_preliminary', 'trialOrder', 'nNeurons', 'nFramesPerTrial',...
        'nTrials', 'nArms', 'A', 'subject')
end
    


    
