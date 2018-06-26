clear; clc; close all

addpath(genpath('~/src/NLW/'))
%

%Logic
doGatherNeuronsTogether= 1


if doGatherNeuronsTogether
f = dir('~/2presults/8arm/')
t = cell(1,8) %all traces all mice all trials
for i = 3:numel(f)
    dirName = fullfile(f(i).folder,f(i).name)
    
    %load neuron trial trace data
    load(fullfile(dirName, 'trialtraces.mat'))
    
    for j = 1:numel(t)
        t{j} = [t{j} ; traces{j}];        
    end
    

end
clearvars A doGatherNeuronsTogether dirName i j 
save('~/2presults/8arm_socialObjectAllNeurons.mat')
else
    load('~/2presults/8arm_socialObjectAllNeurons.mat')
end

%% Define Trial Parts
pre = [0 5]
approach = [ 5 15]
interaction = [15 25]
retreat = [25 35]
post = [35 40]

fTimes = [pre; approach; interaction; retreat; post];
sIdx = round(fTimes*30.98); sIdx(1,1) = 1 %Calculate the indices for each trial epoch
t_m = cell(numel(t), size(sIdx,1))

%Define X Axis
secTick = 5
x = (1:sIdx(end))*1/31;
x = x(1:31*secTick:end);
xticks = 1:sIdx(end);
xticks = xticks(1:31*secTick:end);
set(gca,'XTick', xticks, 'XTickLabel',round(x),'FontSize',4)
xlabel 'time(s)'
ylabel 'Neuron #'

%
tO = cat(3., t{1}, t{3}, t{6});
tS = cat(3,t{2}, t{5}, t{7});
tT = [t{8}];
tB = [t{4}];

%% Individual Neurons and Trials Separated

for i = 1:size(X{1},1)
    for j = numel(X)
        
    shadedErrorBar([],mean(squeeze(tB(i,:,:))'),sem(squeeze(tB(i,:,:))'))
        PlotVerticalLines(sIdx(:,2))

     title(['Neuron ' num2str(i)])
     
    
saveas(h,['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm/Figures/TrialType/Neuron' num2str(i)],'pdf')
    
end


%% Figure with all trial types overlaid

for i = 1:size(tO,1)
    h = figure; h.Units = 'inches'; h.Position = [1 1 3.5 1.5]; hold on;
    s = shadedErrorBar([],mean(squeeze(tO(i,:,:))'),sem(squeeze(tO(i,:,:))')) 
    s.mainLine.Color = 'r';  s.patch.FaceColor = 'r';  s.patch.FaceAlpha = 0.2; s.patch.EdgeAlpha = 0;
    
    
    s = shadedErrorBar([],mean(squeeze(tS(i,:,:))'),sem(squeeze(tS(i,:,:))'))
    s.mainLine.Color = 'b';  s.patch.FaceColor = 'b';  s.patch.FaceAlpha = 0.2;s.patch.EdgeAlpha = 0;

    s = shadedErrorBar([],mean(squeeze(tT(i,:,:))'),sem(squeeze(tT(i,:,:))'));
    s.mainLine.Color = 'c';  s.patch.FaceColor = 'c';  s.patch.FaceAlpha = 0.2;s.patch.EdgeAlpha = 0;

        
    s = shadedErrorBar([],mean(squeeze(tB(i,:,:))'),sem(squeeze(tB(i,:,:))'))
    s.mainLine.Color = 'g';  s.patch.FaceColor = 'g';  s.patch.FaceAlpha = 0.2;s.patch.EdgeAlpha = 0;

    PlotVerticalLines(sIdx(:,2))
    
    set(gca,'XTick',xticks,'xticklabels', round(x),'FontSize',4)

    title(['Neuron ' num2str(i)])
    xlabel('Time (s)')
    
    
saveas(h,['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm/Figures/TrialTypeAll/Neuron' num2str(i)],'pdf')
    
end




%% Novelty Effect
%%
%look at first few trials for each stimulus compared to last few...

for i = 1:size(t{1},1)
    h = figure; h.Units = 'inches'; h.Position = [1 1 8 11]; hold on;
    
    for j = 1:8
        subplot(8,1, j) 
        hold on

        plot(mean(squeeze(t{j}(i,:,1:2)),2),'r')
        plot(mean(squeeze(t{j}(i,:,9:10)),2),'k')

        PlotVerticalLines(sIdx(:,2))   
        set(gca,'XTick',xticks,'xticklabels', round(x),'FontSize',6)
    end
        title(['Neuron ' num2str(i)])
        xlabel('Time (s)')
    
saveas(h,['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm/Figures/Novelty/Neuron' num2str(i)],'pdf')
end










































