function [outputArg1,outputArg2] = grin_2p_ExtractTrialTraces(subject)
%GRIN_2P_EXTRACTTRIALTRACES

init_samz
%check that metadata and neuron data are saved correctly in results
if ~exist(fullfile(resultsPath, subject, 'neuron.mat'))
    fprintf('Need to create neuron file with traces from CNMFE')
else
    'Neuron.mat found!' 
end

if ~exist(fullfile(resultsPath, subject, 'info.mat'))
    fprintf('Need to synchronize metadata')
else
    'Info.mat found!'
end

fprintf('Loading neuron and info files...')
load([resultsPath '/' subject '/neuron.mat'])
load([resultsPath '/' subject '/info.mat'])

%
% downsample frame number
if ~exist('f')
    frames = round(meta.frames/(30/Fs))
end


%sometimes the first frame is 0, get rid of it if so
if frames(1) == 0;
    frames = frames(2:end);
end

% Split Data into Trials
traces_preliminary = [];
for i = 1:size(frames,1)/2
    frame_start = frames(2*i-1); %define start frame
    frame_end = frames(2*i); %define end frame for the trial
    
    if frame_end-frame_start == 155
    traces_preliminary(:,:,i) = C(:,frame_start:frame_end); 
    
    elseif frame_end-frame_start == 154
        traces_preliminary(:,:,i) = C(:,frame_start:frame_end+1); %if off by 1 frame
        
    elseif frame_end - frame_start == 156
         traces_preliminary(:,:,i) = C(:,frame_start:frame_end-1); %if off by 1 frame

    else
        fprintf(num2str(i))
        printf([' A potentional trial is skipped'])           
        end
    
end

%remove 0 matrices from traces
traces = [] ; j = 1;
for i = 1:size(traces_preliminary,3)
    if ~(all(traces_preliminary(:,:,i) == 0))
        traces(:,:,j) = traces_preliminary(:,:,i);
        j = j+1;
    else
        'uh oh'
    end 
end

% Clean up data
clearvars -except A C traces Fs resultsPath subject frames
save(fullfile(resultsPath, subject ,'traces.mat'))

fprintf('Done with trace extraction for each trial')

end

