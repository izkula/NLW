clear; clc
init_samz

fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
f = dir('~/2pdata/20180507/')


%An IMPORTANT note on the trial structure%
%Basically, the frames are unreliable to indicate when a trial starts. This
%is because the roration of the 8arm pushes the trial length over the 30 s
%allowable trial length. Thus,even though the stimulation period begins
%correctly, the object is stil moving towards the mouse during the
%'stationary' period. A reliable way to align trials is to use the
%meta.framesAlined(2*trialnumber) to get teh frame when both the (1) stim e
%nds and (2) retreat begins. I have measured many trial values for several
%mice and trials, and after trial 1, the stim is consistently ~310 frames.
%The stationary lenght is consistently ~212 frames. The trial start to
%trial start is consistently 926ish frames. 

%So, I am going to use the end of the stim period as the algining frame
%number. 

nFramesStim = round(10*30.98)
nFramesBeforeStim = 3*31
nFramesBeforeStimEnd =  nFramesStim + nFramesBeforeStim;
nFramesAfterStimEnd = 10*31;
nFramesStationary = 212
nFramesBeforeStationary = nFramesBeforeStim + (nFramesStim - nFramesStationary);
nFramesTrial = nFramesBeforeStimEnd + nFramesAfterStimEnd


counter = 0;
CT = []; ST = []; MT= []; M = {}; AA = {}; cellCounter = []; CT_S = {};
for i = 3:numel(f)
    foldername = fullfile('~/2pdata/20180507/', f(i).name);
    
    load(fullfile(foldername, 'metadata.mat')) %load metadata
    counter = counter + 1
    CC = []; SS = [];
    for j = 1:3
        
        load(fullfile(foldername, fnames{j}))

        C = fillmissing(C','linear')';
        %S = fillmissing(S,'linear');
        CC = [CC; C];
       %  SS = [SS; S];
        AA{counter,j} = A;
        cellCounter(counter,j) = size(C,1);
      
    end

    %split data into trials
    C_t = nan(size(CC,1),nFramesTrial , size(meta.frames,1)/2);
    %S_t = nan(size(SS,1),nFramesTrial, size(meta.frames,1)/2);
    m_t = nan(1,nFramesTrial, size(meta.frames,1)/2);
    t_width = size(C_t,2)
    
    for j = 1:(size(meta.framesAligned,1)/2)  %for each trial      
        try
         stimEndFrame =  meta.framesAligned(2*j) %define trial marker
     
         C_t(:,:,j) = CC(:,stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);
         %S_t(:,:,j) = SS(:,stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);
         
         m_t(:,:,j) = meta.run(stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);
         
        catch
            fprintf('messed up' )
            i
            j
                
        end      
    end
    CC_S{counter} = CC;
    CT_S{counter} = C_t
    CT = [CT;C_t];
   % ST = [ST; S_t];
    MT = [MT; m_t];
    M{i} = meta; 
    
end

cellCounter = cellCounter

clearvars -except CT CT_S  M  MT AA CC_S cellCounter nFramesBeforeStim nFramesStim nFramesStationary nFramesBeforeStationary nFramesBeforeStimEnd nFramesAfterStimEnd nFramesTrial
stimTrials = repmat([ones(1,8) zeros(1,8)],1,4)
objectTrials = repmat([1 0 1 0 0 1 0 0],1,8)
socialTrials = repmat([0 1 0 0 1 0 1 0],1,8)
blankTrials = repmat([0 0 0 1 0 0 0 0],1,8)
save('~/2presults/8armStim/allNeuronsTrial.mat', '-v7.3')