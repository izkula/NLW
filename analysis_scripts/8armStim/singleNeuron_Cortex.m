
clear; clc
init_samz

fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
f = dir('~/2pdata/20180507/')


load('~/2presults/8armStim/allNeuronsTrial.mat')

approachFrames =  308
stationaryFrames = 312
retreatFrames = 278
ITIFrames = 22

stimTrials = repmat([ones(1,8) zeros(1,8)],1,4)
noStimTrials = abs(stimTrials-1)
objectTrials = repmat([1 0 1 0 0 1 0 0],1,8)
socialTrials = repmat([0 1 0 0 1 0 1 0],1,8)
blankTrials = repmat([0 0 0 1 0 0 0 0],1,8)
toyTrials = repmat([0 0 0 0 0 0 0 1],1,8)


% INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
cTO = CT(:,:,find(objectTrials.*noStimTrials));
cTS = CT(:,:,find(socialTrials.*noStimTrials));
cTB = CT(:,:,find(blankTrials.*noStimTrials));
cTT = CT(:,:,find(toyTrials.*noStimTrials));

cTOS = CT(:,:,find(objectTrials.*stimTrials));
cTSS = CT(:,:,find(socialTrials.*stimTrials));
cTBS = CT(:,:,find(blankTrials.*stimTrials));
cTTS = CT(:,:,find(toyTrials.*stimTrials));

%d = {cTO cTS cTB cTT cTOS cTSS cTBS cTTS}
d = {cTO cTOS cTS cTSS} %cTB cTBS cTT cTTS}
%labels = {'object' 'social' 'blank' 'toy' 'stim-object' 'stim-social' 'stim-blank' 'stim-toy'}
labels = {'object' 'stim-object' 'social' 'stim-social' 'blank' 'stim-blank' 'toy' 'stim-toy'}
%


clear CT ST CTSS CTOS CTBS CTTS CTO CTS CTB CTT
%%
% THIS IS A GOOD FUNCTION USE IT FOR PLOTTING
plot_Single_Neuron(d, fig_dir, approachFrames, stationaryFrames,labels)




















