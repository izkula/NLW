%%%%% Pupil analysis for reaction time arousal test
clear all; close all

addpath(genpath('~/src/OEGAnalyze')); colordef white; 


basePath = '/media/Data/'
dataBasePath = fullfile(basePath, 'data/BpodImageData/')
processedDataBasePath = fullfile(basePath, 'processedData')
resultsBasePath = fullfile(basePath, 'results')
bpodDataPath = '~/Dropbox/Bpod_r0_5/Data';

%% Import data
    %%%% {'mouseName', 'protocolName', 'Date08_YEAR', sessionNumber}, ...
data = { ...
%        {'Dummy Subject', 'MantaTrainToneResponse', 'Jul07_2016', 9}, ... % no pupil, looks good, no reversed frames
%          {'mattlbchatm30', 'MantaTrainToneResponse', 'Jul08_2016', 1}, ... % no pupil, looks good, no reversed frames
%           {'mattlbchatm38', 'MantaTrainToneResponse', 'Jul09_2016', 1}, ... % no pupil, looks good, no reversed frames
%           {'mattlbsstm35', 'MantaTrainToneResponse', 'Jul09_2016', 1}, ... % no pupil, looks good, no reversed frames
%           {'mattlbchatm30', 'MantaTrainToneResponse', 'Jul09_2016', 1}, ... % no pupil, looks good, no reversed frames
%           {'mattlbthm33', 'MantaTrainToneResponse', 'Jul09_2016', 1}, ... % no pupil, looks good, no reversed frames
%              {'mattlbthm59', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbthm60', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbthm61', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbthm62', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbthm63', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbchatm56', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbchatm57', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
%             {'mattlbchatm58', 'MantaTrainStimResponse', 'Oct17_2016', 1}, ...
...
             %{'mattlbchatm51', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
             %{'mattlbchatm52', 'MantaTrainStimResponse', 'Oct18_2016', 1}, ...
             %{'mattlbchatm53', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
             %{'mattlbchatm54', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
...             
             %{'mattlbnpym70', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
             %{'mattlbnpym71', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
             %{'mattlbcartm72', 'MantaTrainStimResponse', 'Oct20_2016', 1}; ...
             %{'mattlbcartm73', 'MantaTrainStimResponse', 'Oct18_2016', 1}; ...
             %{'mattlbthm74', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbthm75', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbthm76', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsertm77', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsertm78', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsertm79', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsstm81', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsstm82', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbsstm83', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbthm84', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbthm85', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbthm86', 'MantaTrainStimResponse', 'Dec12_2016', 1}; ...
             %{'mattlbnpym90', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbnpym91', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbcartm94', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbcartm95', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbcartm96', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbcartm97', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm98', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm99', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm100', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm101', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm102', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm103', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm104', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm105', 'MantaTrainStimResponse', 'Jan28_2017', 1}; ...
             %{'mattlbthm106', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbthm107', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbthm108', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbthm109', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbnpym110', 'MantaTrainStimResponse', 'Feb20_2017', 1}; ...
             %{'mattlbnpym111', 'MantaTrainStimResponse', 'Feb20_2017', 1}; ...
             %{'mattlbnpym112', 'MantaTrainStimResponse', 'Feb20_2017', 1}; ...
             %{'mattlbnpym113', 'MantaTrainStimResponse', 'Feb20_2017', 1}; ...
             %{'mattlbchatm114', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbchatm115', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbchatm116', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             %{'mattlbchatm117', 'MantaTrainStimResponse', 'Feb17_2017', 1}; ...
             {'mattlbcartm134', 'MantaTrainStimResponse', 'Mar22_2017', 1}; ...
             {'mattlbcartm135', 'MantaTrainStimResponse', 'Mar22_2017', 1}; ...
             {'mattlbcartm137', 'MantaTrainStimResponse', 'Mar22_2017', 1}; ...
        };


configPupil = false;
importPupil = true;


whichDatasets = [1:numel(data)];
for k = 1:numel(whichDatasets)
    currData = data{k};
    mouseName = currData{1}; protocolName = currData{2}; date = currData{3}; session = currData{4};
    dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)];
    disp(dataName)
    S = load(MakeBpodBehaviorPath(bpodDataPath, mouseName,  protocolName, date, session));
   
    
    if configPupil
        dataDir = MakeBpodImagePath( dataBasePath, mouseName, protocolName, date, session);
        templateDataName = '';
        trialName = ['Trial', num2str(1, '%05d'), '.mp4']
        [croppedPupilFrames, thresh, eyeMask] = GetPupilFrames(dataDir, trialName, ...
                                                            'templateDataName', templateDataName, 'fromBpod', true);                                             
    end
    

    dataDir = MakeBpodImagePath(dataBasePath, mouseName, protocolName, date, session);
    processedDataDir = strrep(dataDir, 'data', 'processedData');
    templateDataName = processedDataDir;
    if importPupil            
        trackPupil = true;
        ImportPupilBpod(dataDir, 'templateDataName', templateDataName, 'doTrackPupil', trackPupil, 'loadStimChan', true);
        
    end
    close all;
end

%% Load/Analyze traces
doAnalyzeTraces = true;

if doAnalyzeTraces
    for k = 1:numel(whichDatasets)
        currData = data{k};
        mouseName = currData{1}; protocolName = currData{2}; date = currData{3}; session = currData{4};
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)];
        dataDir = MakeBpodImagePath(dataBasePath, mouseName, protocolName, date, session);
        processedDataDir = strrep(dataDir, 'data', 'processedData');


        
        fname = fullfile(processedDataDir, 'pupilTraces.mat');
        PP = load(fname);
        
        fname = fullfile(processedDataDir, 'pupilVids.mat');
        VV = load(fname);
        
        for kk = 1:numel(PP.allMantaTimes)
            %figure, plot(PP.allPupilTimes{kk}, PP.allPupilA{kk}) %%% Note: times should be relative to beginning of trial when Nidaq starts recording
            %hold on; PlotVerticalLines(PP.allStimOnTimes{kk}) %%% Note: times should be relative to beginning of trial when Nidaq starts recording
            %hold on; PlotVerticalLines(PP.allStimOffTimes{kk}) %%% Note: times should be relative to beginning of trial when Nidaq starts recording
        end
%         VideoSlider(VV.allPupilOverlayVid{1})    
    end
end


% SaveVideo( VV.allPupilOverlayVid{1}, fps, savePath, scaleBar, cmap, doColorbar )


