clear; close all; clc

% FOR A SESSIONS WITH 10 STIMS FROM m845 

baseFilePath = '~/2pdata/20180417/';

filePaths = {'m845_stim_004' 'm845_stim_005' 'm845_stim_006' 'm845_stim_007' 'm845_stim_008' ...
    'm860_st5' 'm860_stim_002'  'm860_stim_004' 'm860_stim_006' 'm860_stim_007' ...
    'm861_stim_001' 'm861_stim_002'}


if ~exist('~/2presults/20180417/traces.mat')
a = {} ; c = []; c_raw = []; s = []; frames = {}; cTAvg = []; cTs = {}
    for i = 1:numel(filePaths)
                 
        filePaths{i} %print filepath
        
        f = dir(fullfile(baseFilePath, filePaths{i})); %get all files in the dir
        
        %this is a slapdash way to find the mat file with frames
        for j = 1:numel(f)
            if 100 < f(j).bytes && f(j).bytes < 2002
                load(fullfile(baseFilePath,filePaths{i}, f(j).name))
                break
            end
        end
        
        
        %combine all neuron data for the three planes
        z = {'z0' 'z1' 'z2'}; 
        cOneSession = [];
        for ii = 1:3
            %this is for trial averaging
            try
                load(fullfile(baseFilePath, filePaths{i}, ['AVG_' z{ii} '_neuron.mat']))
                c = [c; C];
                s = [s; S];
                c_raw = [c_raw; C_raw];
                cOneSession = [cOneSession; C]; %this is for trial averaging
            end
        end
            
        %save out the frames
        frames{i} = info.frame ;
        
        %how many seconds do we want for trial averaging?
        nFramesBefore = 4*31;
        nFramesAfter = 4*31;
        nFramesStim = 2*31;
        framesUsed = nFramesBefore + nFramesAfter+nFramesStim;
        
        %Calculate trial averaged calcium traces
            cTOneSession = zeros(size(cOneSession,1), framesUsed, 0.5 * size(info.frame,1));
            for k = 1:size(cTOneSession,3)
                k
                cTOneSession(:,:,k) = ...
                    cOneSession(:,...
                    info.frame(2*k-1)-nFramesBefore : info.frame(2*k-1)+nFramesStim+nFramesAfter-1);
            end
        
            cTs{i} = cTOneSession;
            %get trial average
            cTAvgOneSession = mean(cTOneSession,3)
            cTAvg = [cTAvg ; cTAvgOneSession];
        
   end
        
clearvars -except a c c_raw s frames filepath cTs cTAvg nFramesBefore nFramesAfter nFramesStim framesUsed

save('~/2presults/20180417/traces.mat')





fprintf('resaved')
else

% Load the Data
load('~/2presults/20180417/traces.mat')

end



%% FIGURES
%%%plot showing cortical stim
[~, ind] = max(cTAvg,[],2); %get max index
[~,rsort] = sort(ind); %sort based on this

%sort trial average, calcium and spikes
cTAvgSorted = cTAvg(rsort,:);
cSorted = c(rsort,:);
sSorted = s(rsort,:);

%subtract mean
doMeanSubtract = 0
if doMeanSubtract
    cTAvgSorted = cTAvgSorted - mean(cTAvgSorted(:,1:nFramesBefore),2);
end


doNormalize = 0
if doNormalize
    cTAvgSorted = cTAvgSorted - min(cTAvgSorted, [], 2);
    cTAvgSorted = cTAvgSorted ./ max(cTAvgSorted, [],2);
end

doZ = 1
if doZ
    cTAvgSorted = zscore(cTAvgSorted')';
end


%
%%
h2 = figure; h2.Units = 'inches'; h2.Position = [10 10 3.5 9];
hold on;


imagesc(cTAvgSorted)

secTick = 1
x = (1:length(cTAvgSorted))*1/31;
x = x(1:31*secTick:end);
xticks = 1:length(cTAvgSorted)
xticks = xticks(1:31*secTick:end);
set(gca,'XTick', xticks, 'XTickLabel',round(x),'FontSize',4)
xlabel 'time(s)'
ylabel 'Neuron #'

%light on
rectangle('Position', [nFramesBefore size(cTAvgSorted,1)+1 nFramesStim 50],'FaceColor','r')
PlotVerticalLines([nFramesBefore nFramesBefore+nFramesStim])


export_fig '~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/PSTH_allNeurons.eps'


%% PCA
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(cTAvgSorted);
%%
%FIGURE FOR  EXPLAINED VARIANCE % 
h3 = figure; h3.Units = 'inches'; h3.Position = [5 5 1.5 1.5];
plot(cumsum(EXPLAINED),'ok','MarkerSize',3)
axis([-inf 5 0 100])
set(gca,'XTick',1:5,'YTick',[0 25 50 75 100])
xlabel 'PC'
ylabel 'Variance Explained'

export_fig '~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/PC_Variance.eps'


h4 = figure; h4.Units = 'inches'; h4.Position = [7 7 1.5 1.5]
plot(COEFF(:,1:3))
set(gca,'XTick',xticks,'XTickLabel',round(x),'FontSize',4)
xlabel 'time (s)'
ylabel 'Coeffs'
PlotVerticalLines([nFramesBefore nFramesBefore+nFramesStim])
export_fig '~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/PCs.eps'


%% FIGURE FOR % RESPONSES PER STIMULATION










