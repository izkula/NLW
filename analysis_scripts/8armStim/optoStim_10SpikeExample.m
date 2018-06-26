clear; close all; clc

%% FOR A SESSIONS WITH 10 STIMS FROM m845 

filepath = '~/2pdata/20180417/m845_stim_002'

if ~exist(fullfile(filepath, 'traces.mat'))

%load 2p data
load(fullfile(filepath, 'm845_stim_002_000_002.mat'))

stimStarts = info.frame(1:2:end)

%combine all neuron data for the three planes
z = {'z0' 'z1' 'z2'}; a = {} ; c = []; c_raw = []; s = []; 
for i = 1:3
    
load(fullfile(filepath, ['AVG_' z{i} '_neuron.mat']))

a{i} = A;
c = [c; C];
s = [s; S];
c_raw = [c_raw; C_raw];

end

clearvars -except a c c_raw s info filepath

save(fullfile(filepath, 'traces.mat'))

fprintf('resaved')
end


%% Load the Data
load(fullfile(filepath,'traces.mat'))




%% plot showing cortical stim 
%2 sec before, 2 sec after
nFramesBefore = 2*31
nFramesAfter = 2*31
nFramesStim = 2*31
framesUsed = nFramesBefore + nFramesAfter + nFramesStim

cT = zeros(size(c,1), framesUsed, .5*size(info.frame,1)); %trial averaged
for i = 1:10
    i
    cT(:,:,i) = c(:,info.frame(2*i-1)-nFramesBefore:info.frame(2*i-1)+nFramesAfter+nFramesStim-1);
end

%get trial average
cTAvg = mean(cT,3);
cTAvg = cTAvg - mean(cTAvg(:,1:nFramesBefore)')'; % subtract first second of trialaveraged data

[~, ind] = max(cTAvg,[],2); %get max index
[~,rsort] = sort(ind); %sort based on this

%sort trial average, calcium and spikes
cTAvgSorted = cTAvg(rsort,:)
cSorted = c(rsort,:);
sSorted = s(rsort,:);

%%
%Make some figures
h1 = figure; h1.Units = 'inches'; h1.Position = [10 10 2.5 2.5];


imagesc(cSorted);
colormap


%xaxis
secTick = 10; %seconds for ticks
x = (1:length(cSorted))*1/31;
x = x(1:31*secTick:end);
xticks = 1:length(cSorted);
xticks = xticks(1:31*secTick:end);

secTick = 10; %seconds for ticks
x = (1:length(cSorted))*1/31;
x = x(1:31*secTick:end);
xticks = 1:length(cSorted);
xticks = xticks(1:31*secTick:end);

set(gca,'XTick',xticks,'XTickLabel',round(x),'FontSize',4)
xlabel time(s)
ylabel 'Cortical Neuron Number'
%change amount of data shown
axis([0 190*31 -inf inf])

%add stim rectangles
for i= 1:length(info.frame)/2
    rectangle('Position', [info.frame(2*i-1) -10 , info.frame(2*i)-info.frame(2*i-1), 10],'FaceColor','r','LineWidth',0.001)
end

colormapeditor
export_fig '~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/10Spikesm845002.pdf' 

%%
h2 = figure; h2.Units = 'inches'; h2.Position = [10 10 1.5 2.5];
hold on;


plot(cTAvgSorted','k', 'LineWidth',0.01)

secTick = 1
x = (1:length(cTAvgSorted))*1/31;
x = x(1:31*secTick:end);
xticks = 1:length(cTAvgSorted);
xticks = xticks(1:31*secTick:end);
set(gca,'XTick', xticks, 'XTickLabel',round(x),'FontSize',4)
xlabel 'time(s)'
ylabel 'dF/F'

%light on
rectangle('Position', [62 -15 62 3],'FaceColor','r')

export_fig '~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/PSTH_m845_002.pdf'








%%