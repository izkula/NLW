function SummarizeTracesData2p(dataName, imagePath, bpodTrialDataPath, varargin)

global basePath bpodImagePath resultsPath 

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isGNG', false, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('zplane', [], @isnumeric);
p.addParameter('isNLW', true, @islogical);
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isGNG = p.Results.isGNG;
zplane = p.Results.zplane;
isNLW = p.Results.isNLW;
tDownsampleFactor = p.Results.tDownsampleFactor;

%%% Specify which plots to generate
doSimpleSummary = true;
doGNGSummaryTracePlots = false;
doSummaryBehaviorPlots = false; % THESE SHOULD BE ARGUMENTS


if ~isempty(zplane)
    suffix = ['_z', num2str(zplane), '_'];
else
    suffix = '';
end
    
%%% Load traces into matlab                               
processedDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath);
saveDir = strrep(imagePath, bpodImagePath, fullfile(resultsPath, '2p/')); %% i.e. ~/Dropbox/oeg_results/2p/gad2m11/Image2POlfGoNoGo/Dec16_2015/Session1

%%% Load traces
S = load(fullfile(saveDir, ['traces', suffix, '.mat']));
copyfile(bpodTrialDataPath,fullfile(saveDir, 'behavior.mat'));

Ncells = size(S.traces, 1);
Ntrials = size(S.traces, 3);

%%% Load bpod file
D = load(bpodTrialDataPath);
s = Bpod2Struct(D.SessionData);
try
    frameTimes = s.events{1}.BNC2High; % Just use timestamps from the first trial. 
catch
    t_per_frame = 1/30 * tDownsampleFactor;
    warning(['No frame times found: using t_per_frame = ', num2str(t_per_frame)])
    frameTimes = [0:size(S.fullField, 1)]*t_per_frame; 
end
nZPlanesRecorded = 1;  %%% Record at 30 Hz/number of planes
if isNLW
    try
    if S.metaDatas{1}.volscan
        nZPlanesRecorded = numel(S.metaDatas{1}.otwave_um);
    end
    catch
        warning('Metadatas file not currently included in behavioral file')
    end
end
disp(['nZPlanesRecord: ', num2str(nZPlanesRecorded)]);
frameTimes = frameTimes*nZPlanesRecorded;

Nt = min(size(S.traces,2), numel(frameTimes)); % Number of 2p frames
frameTimes = frameTimes(1:Nt); % We double checked that the 2p outputs *extra* frame TTLs at the *end*.


T = frameTimes;
if strfind(dataName, 'RPE')
    stateChanges = s.events{1}.Tup;
    disp('Using RPE trial types')
    success = GNGSessionSuccess(s);
    hitTrials = s.TrialTypes==1 & success == 1;
    hitTrials = hitTrials(1:Ntrials);
    crTrials = s.TrialTypes==2;
    crTrials = crTrials(1:Ntrials);
    faTrials = s.TrialTypes==3;
    faTrials = faTrials(1:Ntrials);
    xlims = [2, 7];
    ylims = [-0.1 0.3];
elseif strfind(dataName, 'SensoryPrecondTrainNLW')
    disp('Using SensoryPrecondTrainNLW trial types')
    stateChanges = s.events{1}.Tup;
    if ~isempty(s.TrialTypes)
        hitTrials = s.TrialTypes==1;
        hitTrials = hitTrials(1:Ntrials);
        crTrials = s.TrialTypes==2;
        crTrials = crTrials(1:Ntrials);
        faTrials = [];
    else
        hitTrials = 1:Ntrials;
        crTrials = [];
        faTrials = [];
    end
    xlims = [0, 10];
    ylims = [-0.2, 0.5];
elseif strfind(dataName, 'TraceFearJustLick')
    disp('Using TraceFearJustLick trial types')
    events1 = s.events{1};
    stateChanges = [0; events1.GlobalTimer1_End; events1.GlobalTimer2_End; events1.GlobalTimer3_End; events1.GlobalTimer4_End; events1.GlobalTimer5_End]
    hitTrials = 1:Ntrials; %%%% Potentially update this to distinguish between tone, no tone, and shock trials....
    crTrials = [];
    faTrials = [];
    xlims = [0, 60];
    ylims = [-0.5, 1];
elseif strfind(dataName, 'ShapeLightGNG')
    disp('Using ShapeLightGNG trial types')
    events1 = s.events{1};
    stateChanges = [0; events1.Tup(1); events1.Tup(3); events1.Tup(4); events1.Tup(5)]
    hitTrials = 1:Ntrials; %%%% Potentially update this to distinguish between tone, no tone, and shock trials....
    crTrials = [];
    faTrials = [];
    xlims = [0, 6];
    ylims = [-0.5, 1];
elseif strfind(dataName, 'ImageNLWTrainStimGoNoGo')
    disp('Using ImageNLWTrainStimGoNoGo trial types')
    events1 = s.events{1};
    stateChanges = [0; events1.Tup(2); events1.Tup(3); events1.Tup(4)]
    hitTrials = 1:Ntrials; %%%% Potentially update this to distinguish between tone, no tone, and shock trials....
    crTrials = [];
    faTrials = [];
    xlims = [0, 10];
    ylims = [-0.5, 1];
elseif strfind(dataName, 'ImageNLWShapeOlfGNG')
    disp('Using ImageNLWShapeOlfGNG trial types')
    Ntrials = numel(s.TrialTypes);
    success = GNGSessionSuccess(s);
    events1 = s.events{1};
    stateChanges = [0; events1.Tup(2); events1.Tup(3); events1.Tup(4)]
    hitTrials = s.TrialTypes==1 & success == 1;
    hitTrials = hitTrials(1:Ntrials);
    crTrials = s.TrialTypes==2 & success == 1;
    crTrials = crTrials(1:Ntrials);
    faTrials = s.TrialTypes==2 & success == 0;
    faTrials = faTrials(1:Ntrials);
    xlims = [2, 7];
    ylims = [-0.1 0.3];
elseif strfind(dataName, 'ImageNLWTrainOlfGNG')
    disp('Using ImageNLWShapeOlfGNG trial types')
    Ntrials = numel(s.TrialTypes);
    success = GNGSessionSuccess(s);
    events1 = s.events{1};
    stateChanges = [0; events1.Tup(2); events1.Tup(3); events1.Tup(4)]
    hitTrials = s.TrialTypes==1 & success == 1;
    hitTrials = hitTrials(1:Ntrials);
    crTrials = s.TrialTypes==2 & success == 1;
    crTrials = crTrials(1:Ntrials);
    faTrials = s.TrialTypes==2 & success == 0;
    faTrials = faTrials(1:Ntrials);
    xlims = [2, 7];
    ylims = [-0.1 0.3];
elseif strfind(dataName, 'ImageNLWReverseOlfGNG')
    disp('Using ImageNLWReverseOlfGNG trial types')
    Ntrials = numel(s.TrialTypes);
    success = GNGSessionSuccess(s);
    events1 = s.events{1};
    stateChanges = [0; events1.Tup(2); events1.Tup(3); events1.Tup(4)]
    hitTrials = s.TrialTypes==2 & success == 1;
    hitTrials = hitTrials(1:Ntrials);
    crTrials = s.TrialTypes==1 & success == 1;
    crTrials = crTrials(1:Ntrials);
    faTrials = s.TrialTypes==1 & success == 0;
    faTrials = faTrials(1:Ntrials);
    xlims = [2, 7];
    ylims = [-0.1 0.3];
else
   success = GNGSessionSuccess(s);
    hitTrials = s.TrialTypes==1 & success == 1;
    hitTrials = hitTrials(1:Ntrials);
    crTrials = s.TrialTypes==2 & success == 1;
    crTrials = crTrials(1:Ntrials);
    faTrials = s.TrialTypes==2 & success == 0;
    faTrials = faTrials(1:Ntrials);
    xlims = [2, 7];
    ylims = [-0.1 0.3];
end



%%%%%%%%%%% Summary trace plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doSimpleSummary
%%% Save out full field traces
% try
doSaveOutFullField = true;
if doSaveOutFullField
    f=figure('position',[100 100 1024 256]); 
    subplot(1,3,1);
    plot(T, S.dfofFullField(1:Nt,hitTrials),'color',[0.7 0.7 0.7]); hold on;
    plot(T,mean(S.dfofFullField(1:Nt,hitTrials),2),'k'); 
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    title('Hit');
    xlim(xlims); ylim(ylims);
    if ~isempty(crTrials) && any(crTrials>0)
    subplot(1,3,2);
    plot(T, S.dfofFullField(1:Nt,crTrials),'color',[0.7 0.7 0.7]); hold on;
    plot(T,mean(S.dfofFullField(1:Nt,crTrials),2),'k');
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    title('CR');
    xlim(xlims); ylim(ylims);
    end
    if ~isempty(faTrials) && any(faTrials>0)
    subplot(1,3,3);
    plot(T, S.dfofFullField(1:Nt,faTrials),'color',[0.7 0.7 0.7]); hold on;
    plot(T,mean(S.dfofFullField(1:Nt,faTrials),2),'k');
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    title('FA');
    xlim(xlims); ylim(ylims);
    end
    export_fig(f, fullfile(saveDir, 'fullField.pdf'));
end
%%% Save out each cell across trials (don't include neuropil correction
%%% here...will do that for later analyses)
doSaveOutEachTrial = true
if doSaveOutEachTrial
    dfoftraces = zeros(size(S.dfoftraces));
    for i = 1:size(S.dfoftraces, 3)
        smoothVal = 10;
        dfoftraces(:, :, i) = SmoothTraces(S.dfoftraces(:,:,i)', smoothVal)';
    end
    indivTrialDir = fullfile(saveDir, ['trials', suffix]);
    mkdir(indivTrialDir); 
    for i = 1:Ntrials-1
        close all
        f=figure;
        temp = S.dfoftraces(:,:,i);
        temp(isnan(temp)) = 0;
        goodCells = find(std(temp, 0, 2) > 0.03);
        smoothDfof = SmoothRows(S.dfoftraces(goodCells, 1:Nt, i)', 10)';
        smoothDfof = smoothDfof - repmat(mean(smoothDfof, 1), size(smoothDfof,1), 1);
        [~, mind] = max(smoothDfof, [], 2);
        [~, sind] = sort(mind);
        doOrdering = false
        if doOrdering
            smoothDfof = smoothDfof(sind, :);
        end
        
        if ~isnan(temp)
%             maxY = std(temp(:))*1.5;
            maxY = std(smoothDfof(:))*1.5;

            subplot(6,2,[1,3,5, 7, 9]);
%             imagesc(T, 1:numel(goodCells), squeeze(S.dfoftraces(goodCells, 1:Nt, i)),[0 maxY]); colormap(hot);
            imagesc(T, 1:numel(goodCells), squeeze(smoothDfof(:,1:Nt)),[0 maxY]); colormap(hot);
            PlotVerticalLines(stateChanges(2:end-1), 0, size(smoothDfof, 1),'--w');        
            title('Hit');
            xlabel('Time [s]'); ylabel('Cells');
            xlim(xlims);
            hold on; subplot(6,2,[11]);
%             errorshade(T,mean(squeeze(dfoftraces(i, 1:Nt, hitTrials))'),sem(squeeze(dfoftraces(i, 1:Nt, hitTrials))')',[0 1 0],'errorAlpha',0.2); hold on;
%             errorshade(T,mean(smoothDfof, 1),sem(smoothDfof'),[0 1 0],'errorAlpha',0.2); hold on;
            plot(T,smoothDfof');

%             plot(T, mean(squeeze(dfoftraces(i, 1:Nt, hitTrials))'));
            xlim(xlims)
%             ylim(ylims)
            if ~isempty(crTrials)  && any(crTrials>0)
                subplot(6,2,[2, 4, 6, 8, 10]);
                imagesc(T, 1:Ntrials, squeeze(dfoftraces(i, 1:Nt, crTrials))', [0 maxY]); colormap(hot);
                PlotVerticalLines(stateChanges(2:end-1), 0, Ntrials+1,'--w');        
                title('CR');
                xlabel('Time [s]'); ylabel('Trial');  %title(['Cell ', num2str(i)])
                xlim(xlims);
            end
            export_fig(f, fullfile(indivTrialDir, ['cell',num2str(i),'.pdf']));
        end
    end
end


doSaveOutEachCell = false
if doSaveOutEachCell
    dfoftraces = zeros(size(S.dfoftraces));
    for i = 1:size(S.dfoftraces, 3)
        smoothVal = 10;
        dfoftraces(:, :, i) = SmoothTraces(S.dfoftraces(:,:,i)', smoothVal)';
    end
    indivCellDir = fullfile(saveDir, 'cells');
    mkdir(indivCellDir); 
    for i = 1:Ncells
        close all
        f=figure;
        temp = S.dfoftraces(i,:,:);
        if ~isnan(temp)
            maxY = std(temp(:))*1.5;

            subplot(6,2,[1,3,5, 7, 9]);
            imagesc(T, 1:Ntrials, squeeze(dfoftraces(i, 1:Nt, hitTrials))',[0 maxY]); colormap(hot);
            PlotVerticalLines(stateChanges(2:end-1), 0, Ntrials+1,'--w');        
            title('Hit');
            xlabel('Time [s]'); ylabel('Trial');
            xlim(xlims);
            subplot(6,2,[11]);
            errorshade(T,mean(squeeze(dfoftraces(i, 1:Nt, hitTrials))'),sem(squeeze(dfoftraces(i, 1:Nt, hitTrials))')',[0 1 0],'errorAlpha',0.2); hold on;

            plot(T, mean(squeeze(dfoftraces(i, 1:Nt, hitTrials))'));
            xlim(xlims)
            ylim(ylims)
            if ~isempty(crTrials)
                subplot(6,2,[2, 4, 6, 8, 10]);
                imagesc(T, 1:Ntrials, squeeze(dfoftraces(i, 1:Nt, crTrials))', [0 maxY]); colormap(hot);
                PlotVerticalLines(stateChanges(2:end-1), 0, Ntrials+1,'--w');        
                title('CR');
                xlabel('Time [s]'); ylabel('Trial');  %title(['Cell ', num2str(i)])
                xlim(xlims);
            end
            export_fig(f, fullfile(indivCellDir, ['cell',num2str(i),'.pdf']));
        end
    end
end
end

if doGNGSummaryTracePlots
    close all
    %%% Save average cell activity across all cells on each trial
    %avgCellTraces = mean(S.traces, 3);
    %for i = 1:Ncells
        %avgCellTraces(i, :) = smooth(avgCellTraces(i, :),10, 'sgolay');
    %end
    avgCellHit = mean(S.traces(:,:,hitTrials), 3);
    avgCellCR = mean(S.traces(:,:,crTrials),3);
    avgCellFA = mean(S.traces(:,:,faTrials),3);

    figure('position',[100 100 1024 256]);
    subplot(1,3,1);
    plot(repmat(T, Ncells, 1)', avgCellHit','color',[0.7,0.7,0.7]); hold on; %xlim(xlims); 
    plot(T, mean(avgCellHit)','k');
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    ylim([-0.1 0.2]); xlim(xlims);
    title('Hit');
    subplot(1,3,2);
    plot(repmat(T, Ncells, 1)', avgCellCR', 'color', [0.7, 0.7, 0.7]); hold on;
    plot(T, mean(avgCellCR)','k');
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    ylim([-0.1 0.2]); xlim(xlims);
    title('CR');
    subplot(1,3,3);
    plot(repmat(T, Ncells, 1)', avgCellFA', 'color', [0.7, 0.7, 0.7]); hold on;
    plot(T, mean(avgCellFA)','k');
    PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
    ylim([-0.1 0.2]); xlim(xlims);
    title('FA');

    xlabel('Time [s]'); ylabel('df/f');  
    export_fig(f, fullfile(saveDir, ['avgTrialTraces_line.pdf']));


    %%% Now with z-scored traces
    %%% Save average cell activity across all cells on each trial

    disp('Now plotting z-scored')
    for i = 1:Ncells
        for k = 1:size(S.traces, 3)
            avgZCellTraces(i, :, k) = zscore(S.traces(i, :, k));
        end
    end
    avgCellHit = mean(avgZCellTraces(:,:,hitTrials), 3);
    avgCellCR = mean(avgZCellTraces(:,:,crTrials),3);
    avgCellFA = mean(avgZCellTraces(:,:,faTrials),3);

    figure('position',[100 100 1024 256]);
    subplot(1,3,1);
    plot(repmat(T, Ncells, 1)', avgCellHit','color',[0.7,0.7,0.7]); hold on; %xlim(xlims); 
    plot(T, mean(avgCellHit)','k');
    PlotVerticalLines(stateChanges(2:end-2), -1, 1,'--k');
    ylim([-1, 1]); xlim(xlims);
    title('Hit');
    subplot(1,3,2);
    plot(repmat(T, Ncells, 1)', avgCellCR', 'color', [0.7, 0.7, 0.7]); hold on;
    plot(T, mean(avgCellCR)','k');
    PlotVerticalLines(stateChanges(2:end-2), -1, 1,'--k');
    ylim([-1, 1]); xlim(xlims);
    title('CR');
    subplot(1,3,3);
    plot(repmat(T, Ncells, 1)', avgCellFA', 'color', [0.7, 0.7, 0.7]); hold on;
    plot(T, mean(avgCellFA)','k');
    PlotVerticalLines(stateChanges(2:end-2), -1, 1,'--k');
    ylim([-1, 1]); xlim(xlims);
    title('FA');

    xlabel('Time [s]'); ylabel('df/f');  
    export_fig(f, fullfile(saveDir, ['zscore_avgTrialTraces_line.pdf']));

    %%% Now plot zscore heatmaps
    figure('position', [100 100 256 256]);
    imagesc(T, 1:Ncells, avgCellHit);
    xlim(xlims);
    export_fig(f, fullfile(saveDir, ['allcells_zscore_hit.pdf']));

    figure('position', [100 100 256 256]);
    imagesc(T, 1:Ncells,avgCellCR);
    xlim(xlims);
    export_fig(f, fullfile(saveDir, ['allcells_zscore_cr.pdf']));

    figure('position', [100 100 256 256]);
    imagesc(T, 1:Ncells,avgCellFA);
    xlim(xlims);
    export_fig(f, fullfile(saveDir, ['allcells_zscore_fa.pdf']));
    % catch e
    % end
end
%%%%%%%%%%% End summary trace plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Plot bpod trial data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doSummaryBehaviorPlots
    f = figure()
    PlotLicks(s, true);
    PlotVerticalLines(stateChanges(2:end-1), 0, Ntrials+1);
    title('Licks vs trial'), ylabel('Trial'), xlabel('Time [s]')
    export_fig(f, fullfile(saveDir, ['licks.pdf']));

    %%%%%%%%%%% Plot licks and traces overlaid %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%% Plot lick performance %%%%%%%%%%%%%%%%%%%%

    d = strsplit(dataName, '_')
    sessNumber = d{5};
    data = {{d{1}, d{2}, [d{3}, '_', d{4}], sessNumber(end)}};

    [ allTrialData, allFullField, frameTimes, allTrialTypes, ...
        mouseID, allTrialLabels, allLicks, allSuccess, trialNums, ...
        rewardTimes, punishTimes, frameOdorTimes] = GNGLoadAllTrialTraceData2p(fullfile(resultsPath, '2p/'), data);

    figure('position',[100 100 400 150]);
    binSize = 0.1;
    allLickRates = EventRate(allLicks, binSize, 8)';
    if strfind(dataName, 'RPE')
        disp('Using RPE trial types')
        trials1 = allTrialTypes == 1 & allSuccess == 1;
        trials2 = allTrialTypes == 2;
        trials3 = allTrialTypes == 3;

    else
        trials1 = allTrialTypes == 1 & allSuccess == 1;
        trials2 =allTrialTypes == 2 & allSuccess == 1;
        trials3 = allTrialTypes == 2 & allSuccess == 0;
    end
        goLicks = EventRate(allLicks(trials1),binSize,8);
        nogoLicks = EventRate(allLicks(trials2),binSize,8);
        errorNogoLicks = EventRate(allLicks(trials3), binSize,8);

    errorshade((1:size(goLicks,1))*binSize,mean(goLicks,2),sem(goLicks'),[0 1 0],'errorAlpha',0.2); hold on;
    errorshade((1:size(nogoLicks,1))*binSize,mean(nogoLicks,2),sem(nogoLicks'),[0 0 1],'errorAlpha',0.2);
    errorshade((1:size(goLicks,1))*binSize,mean(errorNogoLicks,2),sem(errorNogoLicks'),[1 0 0],'errorAlpha',0.2); hold on;
    legend({'Hit','Correct Reject', 'False Alarm'});
    % PlotVLine(3.3,'k--'); PlotVLine(4.3,'k--');  PlotVLine(4.8,'k--'); PlotVLine(5.8,'k--');
    PlotVLine(stateChanges(2), 'k--'); PlotVLine(stateChanges(3), 'k--'); PlotVLine(stateChanges(4), 'k--');
    % PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k'); hold on;
    xlim([2.3 8]); 
    %ylim([0 12]); set(gca,'FontSize',16);
    xlabel('Time (s)'); ylabel('Lick/s');
    set(gcf,'Color','w');
    saveFigure(fullfile(saveDir, 'avg_licks.pdf'));


    %%%% Plot mean across trials for each cells, with errorshade across cells %%%%%

    figure('position',[100 100 1000 1500]);
    k = 1;
    offset = 3;
    ymin = 1; ymax = 81;
    %T = (ymin:ymax)/15;
    T = frameTimes{1};
    T = T(1:Nt);
    trialNumCutoff=1000;
    currXlim = [-1 3];

    subplot(3,2,1);
    trials = trials1 & trialNums < trialNumCutoff;
    currData = squeeze(mean(S.traces(:,:,trials(1:end-1)), 3)); %%% Mean across trials for each cell
    errorshade(T-offset, squeeze(mean(currData,1)), sem(currData)', GetCmapForGeno(dataName,true), 'errorAlpha', 0.2); hold on;
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    xlim([-1 4]); ylim([-0.1 0.2]);
    licks = EventRate(allLicks(trials),binSize,8);
    Tlicks = (1:size(licks,1))*binSize;
    licks = licks';
    subplot(3,2,2);
    errorshade(Tlicks-offset, mean(licks,1), sem(licks)', [0 0 0], 'errorAlpha', 0.2); hold on;
    title('mean for each cell, error across cells')
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    ylim([0 30]); xlim([-1 4]);
    xlabel('Time [s]'); ylabel('df/f');  
    set(gca,'TickDir','out');

    subplot(3,2,3);
    trials = trials2 & trialNums < trialNumCutoff;
    currData = squeeze(mean(S.traces(:,:,trials(1:end-1)), 3));
    errorshade(T-offset, squeeze(mean(currData,1)), sem(currData)', GetCmapForGeno(dataName,true), 'errorAlpha', 0.2); hold on;
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    xlim([-1 4]); ylim([-0.1 0.2]);
    licks = EventRate(allLicks(trials),binSize,8);
    Tlicks = (1:size(licks,1))*binSize;
    licks = licks';
    subplot(3,2,4)
    errorshade(Tlicks-offset, mean(licks./1000,1), sem(licks./1000)', [0 0 0], 'errorAlpha', 0.2); hold on;
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    ylim([0 30]); xlim([-1 4]);
    xlabel('Time [s]'); ylabel('df/f');  
    set(gca,'TickDir','out');

    subplot(3,2,5);
    trials = trials3 & trialNums < trialNumCutoff;
    currData = squeeze(mean(S.traces(:,:,trials(1:end-1)), 3));
    errorshade(T-offset, squeeze(mean(currData,1)), sem(currData)', GetCmapForGeno(dataName,true), 'errorAlpha', 0.2); hold on;
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    xlim([-1 4]); ylim([-0.1 0.2]);
    licks = EventRate(allLicks(trials),binSize,8);
    Tlicks = (1:size(licks,1))*binSize;
    licks = licks';
    subplot(3,2,6)
    errorshade(Tlicks-offset, mean(licks./1000,1), sem(licks./1000)', [0 0 0], 'errorAlpha', 0.2); hold on;
    PlotVLine(stateChanges(2)-offset, 'k--'); PlotVLine(stateChanges(3)-offset, 'k--'); PlotVLine(stateChanges(4)-offset, 'k--');
    ylim([0 30]); xlim([-1 4]);
    xlabel('Time [s]'); ylabel('df/f');  
    set(gca,'TickDir','out');

    set(gcf,'Color','w');
    saveFigure(fullfile(saveDir, 'meanTraceAndLicks.pdf'));
    end

%%%% To do: Plot full field vs licks %%%%%

%%%% To do: Plot indiv cells vs licks %%%%%



% %%%% Plot 5 highest variance cells, with errorshade across trials %%%%%
% 
% avgCellHit = mean(S.traces(:,:,hitTrials), 3);
% avgCellCR = mean(S.traces(:,:,crTrials),3);
% avgCellFA = mean(S.traces(:,:,faTrials),3);
% 
% maxes = max(avgCellHit, [], 2);
% [mm, mRank] = sort(maxes, 'descend');
% 
% figure('position',[100 100 1000 300]);
% k = 1;
% 
% ymin = 1; ymax = 81;
% %T = (ymin:ymax)/15;
% T = frameTimes{1};
% T = T(1:Nt);
% trialNumCutoff=1000;
% currXlim = [-1 3];
% nTopCells = 5;
% for i=1:nTopCells
%     subplot(3,nTopCells,k);
%     trials = allTrialTypes == 1 & allSuccess == 1 & trialNums < trialNumCutoff;
% 
%     currData = squeeze(S.traces(mRank(i),:,hitTrials));
%     errorshade(T, squeeze(mean(currData,2)), sem(currData'), GetCmapForGeno(dataName,true), 'errorAlpha', 0.2); hold on;
%     PlotVerticalLines(stateChanges(2:end-1), -.1, .3,'--k');
% %     ylim(ylims); xlim(xlims);
% 
% %     currData = squeeze(allTrialData(trials,i,offset:end));
% %     errorshade(T, mean(currData,1), sem(currData)', magentaC, 'errorAlpha', 0.2); hold on;
% 
%     
%     licks = EventRate(allLicks(trials),binSize,8);
%     Tlicks = (1:size(licks,1))*binSize;
%     licks = licks';
%    % licks = licks(Tlicks >= min(T) & Tlicks <= max(T)+1e-2,:)';
%    % Tlicks = Tlicks(Tlicks >= min(T) & Tlicks <= max(T)+1e-2);
% 
%     errorshade(Tlicks, mean(licks/10,1), sem(licks/10)', [0 0 0], 'errorAlpha', 0.2); hold on;
% 
%     trials = allTrialLabels==0 & allTrialTypes == 1 & allSuccess == 1 & trialNums < trialNumCutoff & mouseID == mouse;
%     currData = squeeze(allTrialData(trials,i,offset:end));
%     errorshade(T, mean(currData,1), sem(currData)', greenC, 'errorAlpha', 0.2); hold on;
%     title(L.labels(useRegions(i)));
%     PlotVLine(0,'k--'); PlotVLine(1,'k--'); PlotVLine(1.5,'k--');
%     set(gca,'TickDir','out');
%     xlabel('Time [s]'); ylabel('df/f');  
%     ylim([-0.5 2]);
%     xlim(currXlim);
%     k = k + 1;
% end
end

