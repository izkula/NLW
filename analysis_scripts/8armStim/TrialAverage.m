
%clear; clc; close all
tic
load('~/2presults/8armStim/allNeuronsTrial.mat')
toc
%

fig_dir = '~/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
%%Calculate Z Score Across Entire Data Set
CT_reshaped = reshape(CT,...
    size(CT,1),size(CT,2)*size(CT,3));

CTZ = nanzscore(CT_reshaped')';
CTZ = reshape(CTZ, size(CT,1), size(CT,2), size(CT,3));


%INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
cTO = CTZ(:,:,find(objectTrials.*abs(stimTrials-1)));
cTS = CTZ(:,:,find(socialTrials.*abs(stimTrials-1)));
cTOS = CTZ(:,:,find(objectTrials.*stimTrials));
cTSS = CTZ(:,:,find(socialTrials.*stimTrials));


d = {cTO cTOS cTS cTSS};
labels = {'object' 'Obj stim' 'socal' 'soc stim'}

r = size(CT,1); col = 50;

d_m = [mean(d{1},3) nan(r,col)  mean(d{2},3) nan(r,col) mean(d{3},3) nan(r,col) mean(d{4},3)];


%% Trial Averages 
%Object ON
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

y = d_m;
[m,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');

I(m < 1) = NaN; 

[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 nFramesTrial*2+col 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'ActiveObject'),'pdf')


%% Object OFF
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

y = d_m;
[m,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');

I(m > 0) = NaN; 
[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 nFramesTrial*2+col 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'InActiveObject'),'pdf')


%% Social ON
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 1) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([firstFrame lastFrame])
PlotVerticalLines(nFramesTrial+col+ [firstFrame lastFrame])

axis([2*(nFramesTrial+col) inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'ActiveSocial'),'pdf')

%% Social ODD
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([firstFrame lastFrame])
PlotVerticalLines(nFramesTrial+col+ [firstFrame lastFrame])

axis([2*(nFramesTrial+col) inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'InActiveSocial'),'pdf')

%%  Object  Social ON cells
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 1) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd+ col), y(:,1: nFramesTrial)],[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'ActiveObject_Social'),'pdf')

%% Social Object OFF
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd+ col), y(:,1: nFramesTrial)],[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'AInctiveObject_Social'),'pdf')


%%
%% Social Stim  Object Stim
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 3*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 1) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd), ...
    y(:,  nFramesTrial: 2 * nFramesTrial + col )],[0 2])



PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'ActiveStimObject_Social'),'pdf')

%% Social Stim  Object Stim OFF
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 3*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd), ...
    y(:,  nFramesTrial: 2 * nFramesTrial + col )],[0 2])



PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'AInctiveStimObject_Social'),'pdf')








%% Venn Diagrams
%object stim no stim
y = d_m;
[mNoStim,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim<1) = NaN

[mStim,I] = max(y(:,...
    nFramesTrial+col+nFramesBeforeStim:...
    nFramesTrial+col+nFramesBeforeStim+nFramesStim)');
mStim(mStim<1) = NaN

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenObjectResponsive'],'pdf');

%object stim no stim unresponsive
y = d_m;
[mNoStim,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim>0) = NaN

[mStim,I] = max(y(:,...
    nFramesTrial+col+nFramesBeforeStim:...
    nFramesTrial+col+nFramesBeforeStim+nFramesStim)');
mStim(mStim>0) = NaN

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenObjectNotResponsive'],'pdf');

%% social stim no stim
y = d_m;
[mNoStim,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim<1) = NaN;

[mStim,I] = max(y(:,...
    3*(nFramesTrial+col)+nFramesBeforeStim :3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mStim(mStim<1) = NaN;

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenSocialResponsive'],'pdf')

y = d_m;
[mNoStim,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim>0) = NaN;

[mStim,I] = max(y(:,...
    3*(nFramesTrial+col)+nFramesBeforeStim :3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mStim(mStim>0) = NaN;

respNeuronsNoStim = sum(~isnan(mNoStim));
respNeuronsStim = sum(~isnan(mStim));
respNeuronsBoth = sum(~isnan(mNoStim-mStim));

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenSocialNotResponsive'],'pdf')

%% social object no stim
y = d_m;
[mNoStim,I] = max(y(:,...
        2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mNoStim(mNoStim<1) = NaN;

[mStim,I] = max(y(:,...
    nFramesBeforeStim : nFramesBeforeStim+nFramesStim)');

mStim(mStim<1) = NaN;

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenSocialObjectNoStimResponsive'],'pdf')

%
y = d_m;
[mNoStim,I] = max(y(:,...
        2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mNoStim(mNoStim> 0) = NaN;

[mStim,I] = max(y(:,...
    nFramesBeforeStim : nFramesBeforeStim+nFramesStim)');

mStim(mStim>0) = NaN;

norespNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([norespNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenSocialObjectNoStimNotResponsive'],'pdf')

%% social stim object stim
y = d_m;
[mNoStim,I] = max(y(:,...
    3*(nFramesTrial+col)+nFramesBeforeStim :3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim<1) = NaN;

[mStim,I] = max(y(:,...
    1*(nFramesTrial+col)+nFramesBeforeStim :1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mStim(mStim<1) = NaN;

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenSocialObjectStimResponsive'],'pdf')

y = d_m;
[mNoStim,I] = max(y(:,...
    3*(nFramesTrial+col)+nFramesBeforeStim :3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim>0) = NaN;

[mStim,I] = max(y(:,...
    1*(nFramesTrial+col)+nFramesBeforeStim :1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mStim(mStim>0) = NaN;

notrespNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([notrespNeuronsNoStim respNeuronsStim] , respNeuronsBoth);

saveas(h101, [fig_dir '/VenSocialObjectStimNotResponsive'],'pdf')





