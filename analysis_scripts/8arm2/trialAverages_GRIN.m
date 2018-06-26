clear; clc; close all
tic
load('~/2presults/8arm2/allNeuronsTrial.mat')
toc
%

fig_dir = '~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm2/Figures/'
%%Calculate Z Score Across Entire Data Set
CT_reshaped = reshape(CT,...
    size(CT,1),size(CT,2)*size(CT,3));

CTZ = nanzscore(CT_reshaped')';
CTZ = reshape(CTZ, size(CT,1), size(CT,2), size(CT,3));


%INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
cTO = CTZ(:,:,find(objectTrials));
cTS = CTZ(:,:,find(socialTrials));
cTB = CTZ(:,:,find(blankTrials));
cTT = CTZ(:,:,find(toyTrials));


d = {cTO cTS cTB cTT};
labels = {'object'  'socal ' 'blank' 'toy' }

r = size(CT,1); col = 50;

d_m = [nanmean(d{1},3) nan(r,col)  nanmean(d{2},3) nan(r,col) nanmean(d{3},3) nan(r,col) nanmean(d{4},3)];
%% Trial Averages 
%Object Social
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

y = d_m;
[m,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');

I(m < 0.5) = NaN; 

[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])



axis([1 nFramesTrial*2+col 0 sum(~isnan(I))])
ylabel(' Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'ActiveObject'),'pdf')


%% Object off
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

y = d_m;
[m,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');

I(m > 0) = NaN; 
[~,rsort] = sort(I);
y = y(rsort,:);
imagesc(y,[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])


axis([1 nFramesTrial*2+col 0 sum(~isnan(I))])
ylabel('Unresponsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'InActiveObject'),'pdf')



%%  Social ON cells
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 1*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 0.5) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd+ col), y(:,1: nFramesTrial)],[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'Active_Social_v_Object'),'pdf')

%% Social Object OFF
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 1.5 1.5];

firstFrame = 1*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd+ col), y(:,1: nFramesTrial)],[0 2])
PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])

axis([1 inf 0 sum(~isnan(I))])
ylabel('Trial Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'Inactive_Social_v_Object'),'pdf')


%%
%%  Blank Active and inactive
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 2.25 1.5];

firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 0.5) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim:  lastFrame + nFramesAfterStimEnd + col), ...
    y(:, 1: 2*nFramesTrial + col)],[0 2])


PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([2*nFramesTrial+2*col*1+nFramesBeforeStim 2*nFramesTrial+2*col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])


axis([1 inf 0 sum(~isnan(I))])
ylabel(' Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'Acive _Blank'),'pdf')%%  Blank
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 2.25 1.5];

% % %
firstFrame = 2*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim :  lastFrame + nFramesAfterStimEnd + col), ...
    y(:, 1: 2*nFramesTrial + col)],[0 2])


PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([2*nFramesTrial+2*col*1+nFramesBeforeStim 2*nFramesTrial+2*col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])

axis([1 inf 0 sum(~isnan(I))])
ylabel(' Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'InAcive _Blank'),'pdf')

%% Toy Active and inactive
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 2.25 1.5];

firstFrame = 3*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m < 0.5) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim:  lastFrame + nFramesAfterStimEnd ), nan(322,50) ...
    y(:,  1: 2*nFramesTrial + col)],[0 2])


PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([2*nFramesTrial+2*col*1+nFramesBeforeStim 2*nFramesTrial+2*col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])


axis([1 inf 0 sum(~isnan(I))])
ylabel(' Responsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'Acive _Toy'),'pdf')%%  Blank
h2 = figure; h2.Units = 'Inches'; h2.Position = [1 1 2.25 1.5];

% % %
firstFrame = 3*(nFramesTrial+col)+nFramesBeforeStim
lastFrame =  3*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim

y = d_m;
[m,I] = max(y(:, firstFrame : lastFrame)');

I(m > 0) = NaN; 
I = I
[~,rsort] = sort(I);
y = y(rsort,:);

imagesc([y(:,firstFrame - nFramesBeforeStim:  lastFrame + nFramesAfterStimEnd ), nan(322,50) ...
    y(:,  1: 2*nFramesTrial + col)],[0 2])


PlotVerticalLines([nFramesBeforeStim nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial+col*1+nFramesBeforeStim nFramesTrial+col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([2*nFramesTrial+2*col*1+nFramesBeforeStim 2*nFramesTrial+2*col*1+ nFramesBeforeStim+nFramesStim])
PlotVerticalLines([nFramesTrial + col/2])

axis([1 inf 0 sum(~isnan(I))])
ylabel(' Unresponsive Neuron  #')
xlabel('Time (s)')
set(gca,'XTick',[],'FontSize',6)

saveas(h2,fullfile(fig_dir,'InAcive _Toy'),'pdf')




%% Venn Diagrams
%object stim no stim
y = d_m;
[mNoStim,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim<0.5) = NaN

[mStim,I] = max(y(:,...
    nFramesTrial+col+nFramesBeforeStim:...
    nFramesTrial+col+nFramesBeforeStim+nFramesStim)');
mStim(mStim<0.5) = NaN

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenObjectSocialResponsive'],'pdf');

%object stim no stim unresponsive
y = d_m;
[mNoStim,I] = max(y(:,nFramesBeforeStim:nFramesBeforeStim+nFramesStim)');
mNoStim(mNoStim>0) = NaN;

[mStim,I] = max(y(:,...
    nFramesTrial+col+nFramesBeforeStim:...
    nFramesTrial+col+nFramesBeforeStim+nFramesStim)');
mStim(mStim>0) = NaN;

respNeuronsNoStim = sum(~isnan(mNoStim))
respNeuronsStim = sum(~isnan(mStim))
respNeuronsBoth = sum(~isnan(mNoStim-mStim))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenObjectSocialNotResponsive'],'pdf');

%% blacnk object stim
y = d_m;
[mBlank,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mBlank(mBlank<0.5) = NaN;

[mObj,I] = max(y(:,...
    0*(nFramesTrial+col)+nFramesBeforeStim :0*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mObj(mObj<0.5) = NaN;

respNeuronsNoStim = sum(~isnan(mBlank))
respNeuronsStim = sum(~isnan(mObj))
respNeuronsBoth = sum(~isnan(mBlank-mObj))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenBlankObjectResponsive'],'pdf')
%%
y = d_m;
[mBlank,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mBlank(mBlank<0.5) = NaN;

[mSoc,I] = max(y(:,...
    1*(nFramesTrial+col)+nFramesBeforeStim :1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mSoc(mSoc<0.5) = NaN;

respNeuronsNoStim = sum(~isnan(mBlank))
respNeuronsStim = sum(~isnan(mSoc))
respNeuronsBoth = sum(~isnan(mBlank-mSoc))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

%%
y = d_m;
[mBlank,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mBlank(mBlank>0) = NaN;

[mObj,I] = max(y(:,...
    0*(nFramesTrial+col)+nFramesBeforeStim :0*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mObj(mObj>0) = NaN;

respNeuronsNoStim = sum(~isnan(mBlank))
respNeuronsStim = sum(~isnan(mObj))
respNeuronsBoth = sum(~isnan(mBlank-mObj))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenBlankObjectNotResponsive'],'pdf')
%%
y = d_m;
[mBlank,I] = max(y(:,...
    2*(nFramesTrial+col)+nFramesBeforeStim :2*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');
mBlank(mBlank>0 )= NaN;

[mSoc,I] = max(y(:,...
    1*(nFramesTrial+col)+nFramesBeforeStim :1*(nFramesTrial+col)+ nFramesBeforeStim+nFramesStim)');

mSoc(mSoc>0) = NaN;

respNeuronsNoStim = sum(~isnan(mBlank))
respNeuronsStim = sum(~isnan(mSoc))
respNeuronsBoth = sum(~isnan(mBlank-mSoc))

h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 1.5 0.75];

venn([respNeuronsNoStim respNeuronsStim] , respNeuronsBoth)

saveas(h101, [fig_dir '/VenBlankSocNotResponsive'],'pdf')














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
