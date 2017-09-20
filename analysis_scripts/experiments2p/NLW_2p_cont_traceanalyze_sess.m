%%%% Analyze 2p bpod behavior data, specifcially for volumetric session
%%%% based data (ie not trial based)

clearvars
addpath(genpath('~/git/NLW')); colordef white; 
init_oeg_SAMcomputer;

%%% globals are generated using computer dependent init file (i.e. initOEG or init_oeg_analysis1_NLW)
global basePath bpodDataPath bpodImagePath imagejMacroPath imagejPath resultsPath

%%
%%%%{'mouseName', 'BehaviorProtocol', 'Jan01_2001', SessionNumber}
data = {

%m120 - bReaCheS+
%{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '1', [0,1,2]} %social
%{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '2', [0,1,2]} %object
%{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '3', [0,1,2]} %stim
%{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '4', [0,1,2]} %social2
%{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '5', [0,1,2]} %object
{'m120', 'ImageNLWTrainStimGoNoGo', 'Sep06_2017', '6', [0,1,2]} %stim


};


isNLW = true; %%% Data from Neurolabware microscope
doProvideLoadPath = true;



%% Load deconvolved traces for analysis
shared_results_path = '~/Dropbox/NLW_results/2p/'


for i = 1:numel(data)
    for j = 1:numel(data{i}{5}) %for each plane

        close all
        
        currData = data{i}; %define current data
        
        %name some variables
        mouseName = currData{1}; 
        protocolName = currData{2};  
        currDate = currData{3}; 
        session = currData{4};
        zplane = currData{5}(j); 

    
%dir to save out traces    
saveDir =  MakeBpodImagePath2p(shared_results_path, mouseName, ...
                               protocolName, currDate, session)
%define variable for zstring
zstr = ['z' num2str(zplane)] 

dataType = 'deconv'

%load data for this session and one plane

switch dataType
    case 'traces'
        load(fullfile(saveDir, ['traces_', zstr, '_.mat']))
        d = reshape(permute(traces,[2,1,3]),size(traces,1),[]); %%% the data

        YL = [0,.5] %[0,1]
        outdir = fullfile(saveDir, ['traces_' zstr])
    case 'deconv'
        load(fullfile(saveDir, ['deconvolved_', zstr, '_.mat']), 'dataset')

        d = dataset.img.appendedData.appendedTraces;
        YL = [0,.5]
        outdir = fullfile(saveDir,['deconv_' zstr])
end

    
mkdir(outdir)
t = dataset.img.frameTimes; 
%% Plot individual traces across the whole session
% 
figure, imagesc(t, 1:dataset.Ntrials, sess_across_cells, YL)
xlabel('time [s]')
ylabel('trial #')
xlim([0, max(t(:))])
set(gcf,'color','w');
colorbar()
PlotVerticalLines(reward_times, 1, dataset.Ntrials, 'k-')
title('Mean across cells')
print('-dpdf', fullfile(outdir, 'mean_across_cells.pdf'))




%%% Ah - we weren't properly recording the behavior start times. It looks
%%% as if it is at least 1.5 s after the 2p starts. 
if max(dataset.img.useTrial) == 0
    t = t-2.3; %%% Not 100% sure if this is correct timing, so be aware. 
end
b = dataset.behavior;
event_times = b.events{1}.Tup;
reward_times = b.states{1}.Reward;

%% Average across cells for each trial

mean_across_cells = squeeze(mean(d, 1))';
std_across_cells = squeeze(std(d, 0, 1))';
figure, imagesc(t, 1:dataset.Ntrials, mean_across_cells, YL)
xlabel('time [s]')
ylabel('trial #')
xlim([0, max(t(:))])
set(gcf,'color','w');
colorbar()
PlotVerticalLines(reward_times, 1, dataset.Ntrials, 'k-')
title('Mean across cells')
print('-dpdf', fullfile(outdir, 'mean_across_cells.pdf'))
%% Average across trials for each cell
mean_across_trials = squeeze(mean(d,3));
z_mean_across_trials = squeeze(mean(zscore(d),3));

[~, ind] = max(mean_across_trials(:,20:end), [], 2)
[~, rsort] = sort(ind)
std_across_trials = squeeze(std(d, 0, 3));
figure, imagesc(t, 1:dataset.Ncells, mean_across_trials(rsort, :))
xlabel('time [s]')
% xlim([0, 9])
ylabel('cell #')
set(gcf,'color','w');
colorbar()
PlotVerticalLines(reward_times, 1, dataset.Ncells, 'k-')
% export_fig(fullfile(saveDir, 'mean_across_trials.png'))
title('Mean across trials')
print('-dpdf', fullfile(outdir, 'mean_across_trials.pdf'))
%% Plt individual traces 
doPlotAllTraces = false
if doPlotAllTraces
    figure(203);
    for a=1:size(mean_across_trials,1)
    %     plot(time,dfof(a,:)+(a-1)*5);
        plot(t,mean_across_trials(a,:)+(a-1)*1);
        hold on
    end
    PlotVerticalLines(reward_times)
    xlabel('Time (s');
    ylabel('dF/F');
end
%% Plot overlay traces
figure(204);
for a=1:size(mean_across_trials,1)
%     plot(time,dfof(a,:)+(a-1)*5);
    plot(t,mean_across_trials(a,:));
    hold on
end
ylim(YL)
PlotVerticalLines(reward_times)
xlabel('Time (s');
% xlim([0, 9])
ylabel('dF/F');
title('Average trace for each cell, overlaid')
print('-dpdf', fullfile(outdir, 'overlay_all_traces.pdf'))

%% Plot zscore overlay traces
figure(205);
[~, stim_on_ind] = min(abs(t - reward_times(1)))
for a=1:size(mean_across_trials,1)
%     plot(time,dfof(a,:)+(a-1)*5);
    plot(t,z_mean_across_trials(a,:) - mean(z_mean_across_trials(a,1:stim_on_ind)));
    hold on
end
% ylim(YL)
PlotVerticalLines(reward_times)
xlabel('Time (s');
% xlim([0, 9])
ylabel('dF/F');
title('Average trace for each cell, z-score, baseline-subtracted')
print('-dpdf', fullfile(outdir, 'overlay_all_traces_z.pdf'))
%% Videoslider for each cell
d_cell = permute(d, [3,2,1]);
[~, stim_on_ind] = min(abs(t - reward_times(1)))
[~, stim_off_ind] = min(abs(t - reward_times(2)))


%%% Oh shoot!!! For m594, remember that the timing is not exact because you were not
%%% recording the frame start times. It is delayed...by how much?? 2s?? 
%% Make mask showing locations of 'good cells'
doPlotROIS = true

if doPlotROIS
m = dataset.img.labeledROIs;
good_cells = dataset.goodTraces;

locs = zeros(size(m));
for i = good_cells
    locs(m == i) = i;
end
figure, imagesc(locs)
title('Cell ROIs')
print('-dpdf', fullfile(outdir, 'cell_rois.pdf'))

end

end
end

