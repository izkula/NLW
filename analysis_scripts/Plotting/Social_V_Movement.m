clear; clc; close all

%

%Logic
doGatherNeuronsTogether= 0


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

%%
tO = cat(3., t{1}, t{3}, t{6});
tS = cat(3,t{2}, t{5}, t{7});
tT = [t{8}];
tB = [t{4}];


for i = 1:size(tO,1)
    h = figure; h.Units = 'inches'; h.Position = [1 1 7 10];
    subplot(4,1,1)
    shadedErrorBar([],mean(squeeze(tO(i,:,:))'),sem(squeeze(tO(i,:,:))'))
    PlotVerticalLines(sIdx(:,2))
    
    
    subplot(4,1,2)
    shadedErrorBar([],mean(squeeze(tS(i,:,:))'),sem(squeeze(tS(i,:,:))'))
    PlotVerticalLines(sIdx(:,2))

    subplot(4,1,3)
    shadedErrorBar([],mean(squeeze(tT(i,:,:))'),sem(squeeze(tT(i,:,:))'))
        PlotVerticalLines(sIdx(:,2))

    subplot(4,1,4)
    shadedErrorBar([],mean(squeeze(tB(i,:,:))'),sem(squeeze(tB(i,:,:))'))
        PlotVerticalLines(sIdx(:,2))


     title(['Neuron ' num2str(i)])
     
    
saveas(h,['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm/Figures/TrialType/Neuron' num2str(i)],'pdf')
    
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
for i = 1:numel(t)     
    x = mean(t{i},3);
    
    for j = 1:size(sIdx,1)
    y = mean(x(:,sIdx(j,1):sIdx(j,2)),2);
    t_m{i,j} = y;

    end
end

%%
%Now that we have t_m, we can combine object and social trials together. 
close all
tAAO = [];tAAS = [];
h = []; p = [];
for i = 1:5
tAO = [];
tAS = [];

tAO = [t_m{1,i} t_m{3,i} t_m{6,i}];
tAS = [t_m{2,i} t_m{5,i} t_m{7,i}];

tAAO(:,i) = mean(tAO,2);
tAAS(:,i) = mean(tAS,2);


%ttest
[h(i), p(i)] = ttest(tAAO(:,i), tAAS(:,i));


%barplot
figure;
subplot(2,2,[1 3])
errorb([mean(tAAO(:,i)) mean(tAAS(:,i))], [sem(tAAO(:,i)), sem(tAAS(:,i))])
hold on
bar([mean(tAAO(:,i)) mean(tAAS(:,i))])
%plot([tAAO(:,i) tAAS(:,i)]','k');
hold off
set(gca,'XTick',[1 2],'XTickLabel',{'Object','Social'})
ylabel 'Average Fluorescence'

%scatter plot
subplot(2,2,2)
scatter(tAAO(:,i), tAAS(:,i))
xlabel 'Object'
ylabel 'Social'


%
tAADiff = tAAO(:,i) - tAAS(:,i);
subplot(2,2,4)
histogram(tAADiff)
set(gca,'YScale','log')
ylabel 'n Claustrum Neurons'
xlabel 'Social Preferring .... Object Preferring'


end















%%


close all
%intermediate plots
pos = [1 2 3 4 5 6 10 11]
for i = 1:5
    figure
    for j = [1 3 6 2 5 7 4 8]
        subplot(4,3,pos(j))
        histogram(t_m{j,i})
        set(gca,'YScale', 'log')
    end

end
