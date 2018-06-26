                                                                                                                                                                                        clear; clc
init_samz

fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8arm2/Figures/'
f = dir('~/2pdata/20180522/')


%An IMPORTANT note on the trial structure%
%Basically, the frames are unreliable to indicate when a trial starts. This
%is because the roration of the 8arm pushes the trial length over the 30 s
%allowable trial length. Thus,even though the stimulation period begins
%correctly, the object is stil moving towards the mouse during the
%'stationary' period. A reliable way to align trials is to use the
%meta.framesAlined(2*trialnumber) to get teh frame when both the (1) stim e
%nds and (2) retreat begins. I have measured m any trial values for several
%                                                                    
%mice and trials, and after trial 1, the stim is consistently ~310 frames.
%The stationary lenght is consistently ~212 frames. The trial start to
%trial start is consistently 926ish frames. 

%So, I am going to use the end of the stim period as the algining frame
%number. 

nFramesStationary = 212;
nFramesStim = round(10*30.98);
nFramesBeforeStationary = 3*31
nFramesAfterStationaryEnd = 3*31;
nFramesBeforeStimEnd =   nFramesStim - (nFramesStationary + nFramesBeforeStationary);
nFramesAfterStimEnd = nFramesAfterStationaryEnd
nFramesTrial = nFramesBeforeStationary + nFramesStationary + nFramesAfterStationaryEnd



CT = []; ST = []; MT= []; M = {}; AA = {}; cellCounter = []; CT_S = {}; counter = 1; SSS= {}

for i = 3:numel(f)
    foldername = fullfile('~/2pdata/20180522/', f(i).name);
    if contains(foldername, '8arm2')       
        load(fullfile(foldername, 'metadata.mat')) %load metadata
        CC = []; SS = [];
        for j = 1:3
            try
                load(fullfile(foldername, fnames{j})) %load neurons data

                C = fillmissing(C','linear')';
                S = fillmissing(S,'linear');
                CC = [CC; C];
                SS = [SS; S];
                AA{counter,j} = A;
                cellCounter(counter,j) = size(C,1);
            catch
                fprintf(foldername)
            end
        end

        %split data into trials
        C_t = nan(size(CC,1),nFramesTrial , 64);
        S_t = nan(size(SS,1),nFramesTrial, 64);
        m_t = nan(1,nFramesTrial, 64);
        t_width = size(C_t,2)

        for j = 1:64  %for each trial      
            try
             stimEndFrame =  meta.framesAligned(2*j) %define trial marker

             C_t(:,:,j) = CC(:,stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);
             S_t(:,:,j) = SS(:,stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);

             m_t(:,:,j) = meta.run(stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);

            catch
                fprintf('messed up' )
                i
                j

            end      
        end
        
        SSS{counter} = SS;
        CT_S{counter} = C_t
        CT = [CT;C_t];
        ST = [ST; S_t];
        MT = [MT; m_t];
        M{i} = meta; 
        counter = counter + 1;
    end
end


clearvars -except CT CT_S ST M  MT AA  SSS cellCounter nFramesBeforeStim nFramesStim nFramesBeforeStimEnd nFramesAfterStimEnd nFramesTrial
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8arm2/Figures/'

objectTrials = repmat([1 0 1 0 0 1 0 0],1,8)
socialTrials = repmat([0 1 0 0 1 0 1 0],1,8)
blankTrials = repmat([0 0 0 1 0 0 0 0],1,8)
toyTrials = repmat([0 0 0 0 0 0 0 1],1,8)
save('~/2presults/8arm2/allNeuronsTrial.mat', '-v7.3')

%%





