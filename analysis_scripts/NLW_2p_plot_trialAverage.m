clear; close all; clc
init_samz
%%The script imports trialtraces and then plots all neurons%%

%Logic
doIndividualStimulusFigures = 1

f = dir('~/2presults/8arm/')

%aggregates trace data from all mice and trials and creates one trial
%averaged matrix
trialAvgAll = cell(1,8)
for i = 3:numel(f)
    dirName = fullfile(f(i).folder,f(i).name);
    
    %load neuron trial trace data
    load(fullfile(dirName, 'trialtraces.mat'))
    
    %
    size(traces{1})
    trialAverages = [];
    for ii = 1:numel(traces)
    trialAverages(:,:,ii) = mean(traces{ii},3);
    
    trialAvgAll{ii} = [trialAvgAll{ii}; trialAverages(:,:,ii)];
    end
end


%% group neurons by trial type and calculate normalized value
%%object neurons
objectTrialAvg = mean(cat(3, trialAvgAll{1}, trialAvgAll{3}, trialAvgAll{6}),3);
objectNorm = objectTrialAvg - min(objectTrialAvg, [], 2);
objectNorm = objectNorm ./ max(objectNorm, [], 2);

%group social neurons
socialTrialAvg = mean(cat(3, trialAvgAll{2}, trialAvgAll{5}, trialAvgAll{7}),3);
socialNorm = socialTrialAvg - min(socialTrialAvg, [], 2);
socialNorm = socialNorm ./ max(socialNorm, [], 2);
 
%group empty neurons
emptyTrialAvg = trialAvgAll{4};
emptyNorm = trialAvgAll{4} - min(trialAvgAll{4}, [], 2);
emptyNorm = emptyNorm ./ max(emptyNorm, [], 2);

%group toy mouse neurons
toyTrialAvg = trialAvgAll{8};
toyNorm = trialAvgAll{8} - min(trialAvgAll{8}, [], 2);
toyNorm = toyNorm ./ max(toyNorm, [], 2);

%% reorder by max fluorescence and plot
close all
normTraces = cat(3, objectNorm, socialNorm, emptyNorm, toyNorm);
avgTraces = cat(3, objectTrialAvg, socialTrialAvg, emptyTrialAvg, toyTrialAvg);
trialTypes = {'object', 'social' ,'empty', 'toy'};

doImagePlot = 1;
doLinePlot = 1;


for i = 1:4
d = normTraces(:,:,i);

[~, ind] = max(d,[],2);
[~,rsort] = sort(ind);

d = d(rsort,:);

%time axis
t = linspace(0,size(d,2)/frameSampleRate, size(d,2));

if doImagePlot
    h = figure; h.Units = 'inches'; h.Position = [1 1 5 12]
    imagesc(d)
    title(trialTypes{i})
    ylabel('Neuron Number')
    xlabel('frame(n)')
end


PlotVerticalLines(frameTicks)

%line plot
if doLinePlot
h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 5 12]; hold on;
    for ii = i:size(d,1)
        plot(t,d(ii,:)+ii/5,'r')
    end
    
    PlotVerticalLines(timeTicks)
end


end




















%%
%reorder my max fluorescence in object
[~,ind] = max(objectNorm,[],2);
[~,rsort] = sort(ind);

d = socialNorm(rsort,:);
figure; imagesc(d)
title('Social Ordered by Object Order')
PlotVerticalLines([5*31 5*31+310 5*31+620 length(d)-5*31]);

%reorder by max fluorescence in social
[~,ind] = max(socialNorm,[],2);
[~,rsort] = sort(ind);

d = objectNorm(rsort,:);
figure; imagesc(d)
title('Object Ordered by Social Order')
PlotVerticalLines([5*31 5*31+310 5*31+620 length(d)-5*31]);


%% Put all three types of trials together and then normalize

catTrialAverages = [objectTrialAvg socialTrialAvg emptyTrialAvg];
catTrialNorm = catTrialAverages - min(catTrialAverages, [], 2);
catTrialNorm = catTrialNorm ./ max(catTrialNorm, [], 2);

%reorder by max fluorescence in social
[~,ind] = max(catTrialNorm,[],2);
[~,rsort] = sort(ind);

d = catTrialNorm(rsort,:);
figure; imagesc(d)
title('Normalized across all 3 types')
%PlotVerticalLines([5*31 5*31+310 5*31+620 length(d)-5*31]);

%%


