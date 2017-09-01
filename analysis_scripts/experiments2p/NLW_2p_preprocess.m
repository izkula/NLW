%%%% Preprocess data on the NLW computer before transferring to data
%%%% computer. Only modify this file on the NLW computer. 

% clear all
addpath(genpath('C:\Users\dlab\src\OEGAnalyze')); colordef white; 
init_oeg_NLWcomputer;

%%% globals are generated using computer dependent init file (i.e. initOEG or init_oeg_analysis1_NLW)
global basePath bpodDataPath bpodImagePath imagejMacroPath imagejPath resultsPath

%%
%%%%{'mouseName', 'BehaviorProtocol', 'Jan01_2001', SessionNumber, [zplanes]}
data = {                 
%                {'rbpm12', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 1, [0]};
%                {'rbpm13', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 1, [0, 1, 2]};
%                {'rbpm13', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 2, [0, 1, 2]};
%                {'ckrspm2', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 1, [0]};
%                {'ckrspm2', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 4, [0]};
% % %                {'rbpm14', 'ImageNLWTrainOlfGNG', 'Jul28_2017', 1, [0]};
%                {'rbpm13', 'ImageNLWTrainOlfGNG', 'Jul18_2017', 1, [0, 1, 2]};
%                {'rbpm15', 'ImageNLWTrainOlfGNG', 'Jul31_2017', 1, [0]};
                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Aug05_2017', 2, [0]}
                {'rbpm13', 'ImageNLWTrainOlfGNG', 'Aug06_2017', 5, [0, 1, 2]}
                {'rbpm14', 'ImageNLWTrainOlfGNG', 'Aug06_2017', 2, [0]}
};


isNLW = true; %%% Data from Neurolabware microscope
isContinuousTrial = true;
%%% TO DO TO DO TO DO: Incorporate tDownsampleFactor into registration (so it is
%%% flexible for when doing multiple planes)
tDownsampleFactor = 3; 

%% Merge trials into a single folder if necessary. Do this even if recorded continuously - it puts the data in the right location. 
doAppendTifs = 1

if doAppendTifs
    AppendTifsTogetherDatasets2p(data);
end

%% Motion correct full image stack
doMotionCorrect = 1;

t = cputime;
if doMotionCorrect
    useTemplate = true; %%% Align all trials to the first trial in each session. 
    doDebug = true;
    RegisterMultiTrialDatasets2p(data, 'isContinuousTrial', isContinuousTrial, ...
        'numTrials', [], 'isNLW', isNLW, 'useTemplate', useTemplate, 'doDebug', doDebug, ...
        'tDownsampleFactor', tDownsampleFactor);
end
disp(['Total time for registration:' , num2str(cputime - t)])

%% Generate trial averaged video (no incorporation of trial type)
doMakeAverageVideo = 1;

if doMakeAverageVideo
    MakeAverageVidsDatasets2p(data)
end