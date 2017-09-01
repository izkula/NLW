function MakeAverageVidsData2p(dataName, protocolName, imagePath, behaviorPath, varargin)
% Motion corrects and registers 2p data across trials within a bpod 
% behavioral session.
%
% Saves registered videos to processedData folder. 
%
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath imagejPath imagejMacroPath

%%% Function options
p = inputParser();
p.addParameter('isContinuousTrial', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('zplanes', [], @isnumeric);

p.parse(varargin{:});

isContinuousTrial = p.Results.isContinuousTrial;
numTrials = p.Results.numTrials;
zplanes = p.Results.zplanes;

% processedDataDir = strrep(imagePath, bpodImagePath,  ...
%                           fullfile(basePath, 'processedData', 'BpodImageData'))
processedDataDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath)

if ~exist(processedDataDir, 'dir'); mkdir(processedDataDir); end

%%% Load data
if ~isempty(numTrials)
    subSeq = 1:numTrials
else
    subSeq = []
end


nzplanes = 1;
if ~isempty(zplanes)
    nzplanes = numel(zplanes)
end

B = load(behaviorPath);
b = B.SessionData;
switch protocolName
    case 'ImageNLWReverseOlfGNG'
        go_trial = 2;
        nogo_trial = 1;
        
        ntrials = b.nTrials;
        is_reward = zeros(ntrials, 1);
        is_punish = zeros(ntrials, 1);
        for i = 1:ntrials
            is_reward(i) = ~isnan(b.RawEvents.Trial{i}.States.Reward(1));
            is_punish(i) = ~isnan(b.RawEvents.Trial{i}.States.Punish(1));
        end

        go_trials = (b.TrialTypes == go_trial)';
        nogo_trials = (b.TrialTypes == nogo_trial)';
        hit_trials = (go_trials & is_reward);
        miss_trials = (go_trials & ~is_reward);
        fa_trials = (nogo_trials & is_punish);
        cr_trials = (nogo_trials & ~is_punish);
        success_trials = (hit_trials + cr_trials);
        error_trials = (miss_trials + fa_trials);
        
        trial_lists = [go_trials, nogo_trials, hit_trials, miss_trials, fa_trials, cr_trials, success_trials, error_trials];
        trial_labels = {'go', 'nogo', 'hit', 'miss', 'fa', 'cr', 'success', 'error'}
        
    case 'ImageNLWTrainOlfGNG'
        go_trial = 1;
        nogo_trial = 2;
        
        ntrials = b.nTrials;
        is_reward = zeros(ntrials, 1);
        is_punish = zeros(ntrials, 1);
        for i = 1:ntrials
            is_reward(i) = ~isnan(b.RawEvents.Trial{i}.States.Reward(1));
            is_punish(i) = ~isnan(b.RawEvents.Trial{i}.States.Punish(1));
        end

        go_trials = (b.TrialTypes == go_trial)';
        nogo_trials = (b.TrialTypes == nogo_trial)';
        hit_trials = (go_trials & is_reward);
        miss_trials = (go_trials & ~is_reward);
        fa_trials = (nogo_trials & is_punish);
        cr_trials = (nogo_trials & ~is_punish);
        success_trials = (hit_trials + cr_trials);
        error_trials = (miss_trials + fa_trials);
        
        trial_lists = [go_trials, nogo_trials, hit_trials, miss_trials, fa_trials, cr_trials, success_trials, error_trials];
        trial_labels = {'go', 'nogo', 'hit', 'miss', 'fa', 'cr', 'success', 'error'}
        
    case 'ImageNLWTrainStimGoNoGo'
        ntrials = b.nTrials;
        trial_lists = ones(ntrials, 1);
        trial_labels = {'stim'};
        
        
end
        
    
for ttype = 1:numel(trial_labels)
    for zz = 1:nzplanes

        if ~isempty(zplanes)
            zplane = zplanes(zz);
        else
            zplane = [];
        end

        viddir = GetProcessedDataDir2p(fullfile(imagePath, dataName), basePath, bpodImagePath)
        vidfile = fullfile(['z', num2str(zplane)], 'reg', 'registered_noDenoise.tif')
        outdir = fullfile(viddir, trial_labels{ttype})
        whichTrials = find(trial_lists(:, ttype)); %%%% Can eventually fill this in according to behavior data.

        MakeAverageVid(viddir, vidfile, outdir, whichTrials)
    end
end

end

