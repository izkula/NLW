%%%% Analyze 2p bpod behavior data for *continuously* recorded session (all trials are together in one movie).

clear all

addpath(genpath('~/git/NLW/'));
% addpath(genpath('C:\Users\dlab\src\OEGAnalyze')); 
colordef white; 
% init_oeg_analysis1_NLW;
% init_oeg_NLWcomputer;
init_oeg_SAMcomputer;


%%% globals are generated using computer dependent init file (i.e. initOEG or init_oeg_analysis1_NLW)
global basePath bpodDataPath bpodImagePath imagejMacroPath imagejPath resultsPath

%%
%%%%{'mouseName', 'BehaviorProtocol', 'Jan01_2001', SessionNumber, [zplanes]}
data = {
%     {'ck2RSPm1', 'SensoryPrecondTrainNLW', 'Jan16_2017', 1, []}; ...  %%% Just sound
%     {'ck2RSPm1', 'SensoryPrecondTrainNLW', 'Jan16_2017', 3, []}; ...  %%% Just light
%     {'ck2RSPm1', 'SensoryPrecondTrainNLW', 'Jan16_2017', 4, []}; ...  %%% Just odor
%     {'ckRSPm2', 'SensoryPrecondTrainNLW', 'Jan16_2017', 1, []}; ...  %%% Just sound
%     {'ckRSPm2', 'SensoryPrecondTrainNLW', 'Jan16_2017', 2, []}; ...  %%% Just light
%     {'ckRSPm2', 'SensoryPrecondTrainNLW', 'Jan16_2017', 3, []}; ...  %%% Just odor
%     {'ckrspm1', 'ShapeLightGNG', 'Feb12_2017', 1, [0]}; ... %%% Just flashed light, no reward
%     {'ckrspm2', 'ShapeLightGNG', 'Feb12_2017', 1, [0]}; ... %%% Just flashed light, no reward
%     {'ckRSPm2', 'ShapeLightGNG', 'Feb12_2017', 1, [0]}; ... %%% Just flashed light, no reward
%     {'ckrspm1', 'ShapeLightGNG', 'Feb22_2017', 3, [0]}; ... %%% Wasn't licking 
%     {'ckrspm1', 'ShapeLightGNG', 'Feb22_2017', 5, [0]}; ... %%% Not very thirsty %%% select cells from here on. 
%     %%{'ckRSPm2', 'ShapeLightGNG', 'Feb22_2017', 1, [0]}; ... %%% no licking NEED TO IMPORT
%     {'ckRSPm2', 'ShapeLightGNG', 'Feb22_2017', 4, [0]}; ... %%% no licking
%     {'ckRSPm2', 'ShapeLightGNG', 'Feb22_2017', 5, [0]}; ... %%% licked

% % % %       {'rbpm11', 'ImageNLWTrainLickGoNoGo', 'Apr20_2017', 3, [0]}; ... %%% licked %%% TO ANALYZE
% % % %       {'rbpm12', 'ImageNLWTrainLickGoNoGo', 'Apr21_2017', 3, [0]}; ... %%% licked %%% TO ANALYZE
      
%%%%%%%% CLA stim with RSP image, round 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         {'m594', 'ImageNLWTrainStimGoNoGo', 'May04_2017', 1, [0]};
%         {'m594', 'ImageNLWTrainStimGoNoGo', 'May04_2017', 2, [0]};
%         {'m592', 'ImageNLWTrainStimGoNoGo', 'May04_2017', 1, [0]}; ...
%          {'m592', 'ImageNLWTrainStimGoNoGo', 'May04_2017', 3, [0]}; ...
%         {'m592', 'ImageNLWTrainStimGoNoGo', 'May04_2017', 5, [0]}; ...
%         {'ckrspm2', 'ImageNLWTrainStimGoNoGo', 'Jun01_2017', 1, [0]};

%         {'rbpm12', 'ImageNLWShapeOlfGNG', 'May18_2017', 1, [0]}; ... %%% 
%         {'rbpm12', 'ImageNLWTrainOlfGNG', 'May18_2017', 1, [0]}; ... %%% 


%%%%%%%% CLA stim with RSP image, round 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Note: remember that whisker movement may have led to a lot of the
%%%% activity?
%          {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 2, [0]}; %% 3s on, 3s ITI
%          {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 3, [0]}; %% 3s on, 3s ITI
%          {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 4, [0]}; %% 3s on, 12s ITI
%          {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 5, [0]}; %% 60s on, 60s ITI
%           {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 6, [0]}; %%  3s on, 3s ITI (left)
%           {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 7, [0]}; %%  3s on, 12 s ITI 
%           {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul05_2017', 9, [0]}; %%  60s on, 60s ITI 
% 
%             {'m295',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 1, [0]}; % right 3s on 3s iti
%             {'m295',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 2, [0]}; % right 60s on 60s iti
%             {'m295',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 4, [0]}; % left 3s on 3s iti
%             {'m295',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 5, [0]}; % left 60s on 60s iti
% 
%             {'m118',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 2, [0]}; % right 3s on 3s iti
% % %             {'m118',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 3, [0]}; % right 60s on 60s iti
% 
%             {'m407',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 2, [0]}; % left 3s on 3s iti
%             {'m407',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 3, [0]}; % left 60s on 60s iti
%             {'m407',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 4, [0]}; % right 3s on 3s iti
%             {'m407',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 5, [0]}; % right 60s on 60s iti

%             {'m293',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 1, [0]}; % right 3s on 3s iti
%             {'m293',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 2, [0]}; % right 60s on 6s iti
%             {'m293',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 3, [0]}; % left 3s on 3s iti
%             %%%{'m293',  'ImageNLWTrainStimGoNoGo', 'Jul07_2017', 4, [0]}; % left 60s on 6s iti (doesn't appear to have actually been recorded?)

% 
%                {'rbpm12', 'ImageNLWReverseOlfGNG', 'Jul15_2017', 2, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul15_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul14_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul12_2017', 1, [0]};
%                {'ckrspm2', 'ImageNLWTrainOlfGNG', 'Jul14_2017', 2, [0]};


%             {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 1, [0]};
%             {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 2, [0]};
%             {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 3, [0]};
%             {'m119', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 4, [0]};
%             {'m407', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 1, [0]};
%             {'m407', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 2, [0]};
%             {'m118', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 1, [0]};
%             {'m118', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 2, [0]};
%             {'m295', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 1, [0]};

%             {'m295', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 2, [0]};
%             {'m120', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 1, [0]};
%             {'m120', 'ImageNLWTrainStimGoNoGo', 'Jul13_2017', 2, [0]};

%                {'rbpm12', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 1, [0]};
%                {'rbpm13', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 1, [0, 1, 2]};
%                {'rbpm13', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 1, [0, 1, 2]};
%                {'ckrspm2', 'ImageNLWTrainOlfGNG', 'Jul22_2017', 1, [0]};
%                {'ckrspm2', 'ImageNLWReverseOlfGNG', 'Jul22_2017', 4, [0]};

              % {'rbpm12', 'ImageNLWReverseOlfGNG', 'Jul15_2017', 2, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul15_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul14_2017', 1, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'Jul12_2017', 1, [0]};
%                {'ckrspm2', 'ImageNLWTrainOlfGNG', 'Jul14_2017', 2, [0]};
%                {'rbpm12', 'ImageNLWTrainOlfGNG', 'May18_2017', 1, [0]};

%                 {'rbpm14', 'ImageNLWTrainOlfGNG', 'Jul28_2017', 1, [0]};
%              {'rbpm12', 'ImageNLWTrainOlfGNG', 'Aug05_2017', 2, [0]}
%             {'rbpm13', 'ImageNLWTrainOlfGNG', 'Aug06_2017', 5, [0, 1, 2]}
%             {'rbpm14', 'ImageNLWTrainOlfGNG', 'Aug06_2017', 2, [0]}
            
%              {'rbpm13', 'ImageNLWTrainOlfGNG', 'Jul18_2017', 1, [0, 1, 2]} %%% first day, first round --- z planes [0, 1, 2]

%%%%%%%%%%%%%%%   SAM 2p CLAUSTURUM DATA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%m118 - bReaChes+
{'m118', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '2', [0,1,2]} %social
{'m118', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '3', [0,1,2]} %social2
{'m118', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '4', [0,1,2]} %objects
{'m118', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '5', [0,1,2]} %optostim

%m119 - bReaCheS+
%{'m119',  'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '1', [0,1,2]} 
%{'m119', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '3 ', [0]} 
%{'m119', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '4', [0,1,2]} %social

%m120 - bReaCheS+
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '1', [0,1,2]} %social
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '2', [0,1,2]} %object
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '3', [0,1,2]} %stim
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '4', [0,1,2]} %social2
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '5', [0,1,2]} %object
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '6', [0,1,2]} %stim

%m293
 {'m293', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '1', [0,1,2]} %social
 {'m293', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '2', [0,1,2]} %object
 {'m293', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '3', [0,1,2]} %stim
% 
% %m295
 {'m295', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '1', [0,1,2]} %social
 {'m295', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '2', [0,1,2]} %social2
 {'m295', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '3', [0,1,2]} %object1
 {'m295', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '4', [0,1,2]} %stim
% 
% %m407
 {'m407', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '1', [0,1,2]} %social
 {'m407', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '3', [0,1,2]} %object1
 {'m407', 'ImageNLWTrainStimGoNoGo', 'Sep07_2017', '4', [0,1,2]} %stim
};


isNLW = true; %%% Data from Neurolabware microscope
isContinuousTrial = true;
%%% TO DO TO DO TO DO: Incorporate tDownsampleFactor into registration (so it is
%%% flexible for when doing multiple planes)
tDownsampleFactor = 1; 

%% Merge trials into a single folder if necessary. Do this even if recorded continuously - it puts the data in the right location. 
doAppendTifs = 0
if doAppendTifs
    AppendTifsTogetherDatasets2p(data);
end

%% Motion correct full image stack
doMotionCorrect = 0

t = cputime;
if doMotionCorrect
    useTemplate = true; %%% Align all trials to the first trial in each session. 
    doDebug = true;
    RegisterMultiTrialDatasets2p(data, 'isContinuousTrial', isContinuousTrial, ...
        'numTrials', [], 'isNLW', isNLW, 'useTemplate', useTemplate, 'doDebug', doDebug);
end
disp(['Total time for registration:' , num2str(cputime - t)])

%% Generate trial averaged video (no incorporation of trial type)
doMakeAverageVideo = 0

if doMakeAverageVideo
    MakeAverageVidsDatasets2p(data)
end

%this is the end of NLW_2p_preprocess%

%% Do ICA-based automated cell selection (requires user input)
doICAselection = 1
if doICAselection
    ICAselectionMultiTrialDatasets2p(data)
end

%% Now find cell masks (automatic)
doICAsegmentation = 1

if doICAsegmentation
    doPlotting = 0;
    ICAsegmentationMultiTrialDatasets2p(data, 'doPlotting', doPlotting)
end


%% Now take those cell masks and use existing code to clean it up. (user input)

doAdjustSegmentation = 0

if doAdjustSegmentation 
    defaultCellRadius = 3; % For manually selected neurons (as opposed to automated)
    SelectCellsICAMultiTrialDatasets2p(data, 'cellRadius', defaultCellRadius);
end


%% Extract traces (including neuropil etc). 

doSaveTraces = 1;

if doSaveTraces
    SaveTracesMultiTrialDatasets2p(data, 'numTrials', [], 'isNLW', isNLW, ...
                                    'isContinuousTrial', isContinuousTrial, ...
                                    'tDownsampleFactor', tDownsampleFactor);
end

%%% Load bpod data and summarize traces (takes a while) (use deconvolve traces instead of this)
doSummarizeTraces = 0;

if doSummarizeTraces
    SummarizeTracesDatasets2p(data, 'isNLW', isNLW, 'tDownsampleFactor', tDownsampleFactor);
end

%%% deconvolve traces %%% TO DO TO DO

doDeconvolveTraces =1;
if doDeconvolveTraces
    DeconvolveTracesDatasets(data, 'numTrials', [], 'isNLW', isNLW, 'tDownsampleFactor', tDownsampleFactor);
end

fprintf('Done with processing, check your Dropbox dir')


% % %% Make some plots of the trace data
% % doPlots = false;
% % 
% % if doPlots
% % AA = load('/home/izkula/Dropbox/NLW_results/2p/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session1/traces_z0_.mat')
% % BB = load('/home/izkula/Dropbox/NLW_results/2p/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session1/deconvolved_z0_.mat')
% % dfof = AA.dfoftraces;
% % startFrames = ceil(AA.startFrames/tDownsampleFactor);
% % figure(203);
% % for a=1:size(dfof,1)
% % %     plot(time,dfof(a,:)+(a-1)*5);
% %     plot(dfof(a,:)+(a-1)*5);
% %     hold on
% % end
% % PlotVerticalLines(startFrames)
% % xlabel('Time (s');
% % ylabel('dF/F');
% % 
% % deconv = BB.dataset.img.deconv;
% % startFrames = ceil(BB.dataset.startFrames/tDownsampleFactor);
% % figure(203);
% % for a=1:size(dfof,1)
% % %     plot(time,dfof(a,:)+(a-1)*5);
% %     plot(dfof(a,:)+(a-1)*5);
% %     hold on
% % end
% % PlotVerticalLines(startFrames)
% % xlabel('Time (s');
% % ylabel('dF/F');
% % end


%% Do ICA cell selection
% pathname = '/home/izkula/Data/processedData/BpodImageData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session1/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session1/z0/'
% filedate = 'reg_meansub/'
% filename = 'registered_noDenoise_meansub'
% flims = ''
% nPCs = 200
% dsamp = 2
% badframes = []
% mu = 0.1
% run_pcaica(pathname,filedate,filename,flims, nPCs,dsamp,badframes,mu);

% %% Find cell masks
% smwidth = 0.2
% thresh = 4
% arealims = 20
% plotting = 0
% run_findcellmasks(pathname,filedate,filename,smwidth, thresh, arealims, plotting);
% 
% %% Delete cell masks (need to work on some more)
% run_deletecellmasks(pathname,filedate,filename);


% sum = 
%  figure, imagesc(squeeze(sum(dfof, 1)))

%%% Make trial matrix
% trials = ReshapeToTrials(dfof, startFrames);
% trials = zeros(size(dfof, 1), ceil(max(diff(startFrames))), numel(startFrames));
% for i = 1:numel(startFrames)-1
%     trials(:,1:startFrames(i+1)-startFrames(i), i) = dfof(:, startFrames(i):startFrames(i+1)-1);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %% Plot cell masks (need to work on some more)
% run_plotcellmasks(pathname, filedate,filename);
% 
% %% Extract traces
% frametime = 1/10;
% deconvtauoff = 0.1;
% slidingwind = 10;
% percentfilter = 8;
% lowpassfilt = 3;
% run_applycellmasks(pathname,filedate,filename,frametime,...
%     deconvtauoff,slidingwind,percentfilter,...
%     lowpassfilt);
% 
% run_plotdfof(frametime,pathname,filedate,filename);


%%
% %% Motion correct and register all videos (takes a while)
% 
% doMotionCorrect = false;
%     
% t = cputime;
% if doMotionCorrect
%     useTemplate = true; %%% Align all trials to the first trial in each session. 
%     doDebug = false;
%     RegisterMultiTrialDatasets2p(data, 'numTrials', [], 'isNLW', isNLW, 'useTemplate', useTemplate, 'doDebug', doDebug);
% end
% disp(['Total time for registration:' , num2str(cputime - t)])
% 
% 
% % Save out the stats (i.e. max, skewness) across trials (takes a while)
% doSaveRegisteredStats = false
% if doSaveRegisteredStats
%     SaveRegisteredStatsDatasets2p(data);
% end


% %% Select cells for each session (requires user input, though mostly automated)
% %%% TO DO: SOFTWARE TO SELECT SAME CELLS ACROSS DAYS. Improve cell
% %%% selection!
% doSelectCells = false;   
% 
% if doSelectCells 
%     defaultCellRadius = 6; % For manually selected neurons (as opposed to automated)
%     SelectCellsMultiTrialDatasets2p(data, 'cellRadius', defaultCellRadius);
% end
%     
% 
% doSelectCellsICA = true; %%% Use the Mukamel et al instead. 
% if doSelectCellsICA
%     doAppendTiffs = true
%     SelectCellsICAMultiTrialDatasets2p(data, 'doAppendTiffs', doAppendTiffs);
% end
% 
% 
% %% Save out traces all trials in each session, summarize data (takes a while)
% 
% doSaveTraces = false;
% 
% if doSaveTraces
%     SaveTracesMultiTrialDatasets2p(data, 'numTrials', [], 'isNLW', isNLW);
% end
% 
% %% Load bpod data and summarize traces (takes a while) (use deconvolve traces instead of this)
% 
% doSummarizeTraces = false;
% 
% if doSummarizeTraces
%     SummarizeTracesDatasets2p(data, 'isNLW', isNLW);
% end
% 
% %% deconvolve traces
% doDeconvolveTraces = false;
% if doDeconvolveTraces
%     DeconvolveTracesDatasets(data, 'numTrials', [], 'isNLW', isNLW);
% end
