clear; clc
init_samz

fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/'
f = dir('~/2pdata/20180507/')

CC = []; SS = []; CT = []; ST = []; M = {};
for i = 3:numel(f)
    foldername = fullfile('~/2pdata/20180507/', f(i).name);
    
    load(fullfile(foldername, 'metadata.mat'))
    counter = 1;
    for j = 1:3
        
        load(fullfile(foldername, fnames{j}))

        C = fillmissing(C','linear')';
        S = fillmissing(S,'linear');
        CC = [CC; C];
        SS = [SS; S];
        
        counter = counter + 1;
    end

    %split data into trials
    C_t = zeros(size(CC,8),1082 , size(meta.frames,1)/2);
    S_t = zeros(size(SS,1),1082, size(meta.frames,1)/2);
    t_width = size(C_t,2)
    
    for j = 1:(size(meta.frames,1)/2)        
        try
        C_t(:,:,j) = CC(:,meta.frames(2*j-1) - round(nFramesBaseline) + meta.frameStart : (meta.frames(2*j-1) + + meta.frameStart + 927-1));
        S_t(:,:,j) = SS(:,meta.frames(2*j-1) - round(nFramesBaseline) + meta.frameStart : (meta.frames(2*j-1) + meta.frameStart  + 927-1));
        
        catch %for sessions that are longer than included..
            try
                C_t(:,:,j) = [CC(:,meta.frames(2*j-1) - round(nFramesBaseline) + meta.frameStart : length(CC)) nan(size(CC,1), (meta.frames(2*j-1) + meta.frameStart + 927-1) - 60000)];
                S_t(:,:,j) = [SS(:,meta.frames(2*j-1) - round(nFramesBaseline) + meta.frameStart : length(SS)) nan(size(SS,1), (meta.frames(2*j-1) + meta.frameStart + 927-1) - 60000)];
            catch
                C_t(:,:,j) =  nan(size(CC,1), 1082);
                S_t(:,:,j) =  nan(size(SS,1), 1082);
            end
            
     
        end        
    end
    
    CT = [CT;C_t];
    ST = [ST; S_t];
    M{i} = meta;
    
end

clearvars -except CT ST M

save('~/2presults/8armStim/allNeuronsTrial.mat', '-v7.3')
%% Housekeeping

init_samz
fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/'
load('~/2presults/8armStim/allNeuronsTrial.mat')

stimStartFrame = 9*31;
stimEndFrame = stimStartFrame + 10*31;
stim = repmat([ones(1,8) zeros(1,8)],1,4)
objectTrials = repmat([1 0 1 0 0 1 0 0],1,8)
socialTrials = repmat([0 1 0 0 1 0 1 0],1,8)
blankTrials = repmat([0 0 0 1 0 0 0 0],1,8)

%get one value for each trial from each neuron
CTA = []; STA = []; 
for i = 1:64
    CTA(:,i) = nanmean(CT(:, stimStartFrame:stimEndFrame,i),2);
    STS(:,i) = nansum(ST(:, stimStartFrame:stimEndFrame,i),2);
end
%% INDIVIDUAL NEURON ANLAYSIS -- STIM NO STIM -- USING CALCIUM
%how does the a neuron's activity change with stimulation

% % %for OBJECT trials
COS = (CTA .* stim) .* objectTrials;
COE = (CTA .*abs(stim -1)) .* objectTrials;

COS(:, all(~any(COS),1)) = [];
COE(:, all(~any(COE),1)) = [];

%calcuate a p value for stim vs nostim object trials
p = ones(size(COS,1),1);
for i = 1:size(COS,1)
    p(i) = ranksum(COS(i,:), COE(i,:));
end

%find the cells that are significantly different
I_obj = find(p < 0.01)

% % % for SOCIAL trials
CSS = (CTA .* stim) .* socialTrials;
CSE = (CTA .*abs(stim -1)) .* socialTrials;

CSS(:, all(~any(CSS),1)) = [];
CSE(:, all(~any(CSE),1)) = [];

%calcuate a p value for stim vs nostim object trials
p = ones(size(CSS,1),1);
for i = 1:size(COS,1)
    p(i) = ranksum(CSS(i,:), CSE(i,:));
end

%find the cells that are significantly different
I_soc = find(p < 0.01)




%% Let's plot all of these cells that are different.
plotObjectCells = 1
if plotObjectCells
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 8 2]
hold on
for i = 1:numel(I_obj)
    scatter((i)*ones(1,12), COS(I_obj(i),:),'r','filled')
    scatter((i+0.25)*ones(1,12),  COE(I_obj(i),:),'k','filled')
end

xlabel 'Neuron Numbr'
ylabel 'Average Fl'
end


plotSocialCells = 0
if plotSocialCells
h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 8 2]
hold on
for i = 1:numel(I_soc)
    scatter((i)*ones(1,12), CSS(I_obj(i),:),'r','filled')
    scatter((i+0.25)*ones(1,12),  CSE(I_obj(i),:),'k','filled')
end

xlabel 'Neuron Numbr'
ylabel 'Average Fl'
end

%plot a venn diagram of significantly modulated neurons
plotVennDiagram = 1
if plotVennDiagram
    I_intersect = intersect(I_obj, I_soc);
    h3 = figure; h3.Units = 'inches'; h3.Position = [5 5 2 1]
    venn([length(I_obj) length(I_soc)], length(I_intersect))
    saveas(h3,fullfile(fig_dir, 'Venn'),'pdf')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE ZSCORES--  AND PERFORM UNSUPERVISED LEARNING %%%%
%we only want to include social and object trials
zTrials = [find(objectTrials) find(socialTrials)];

%covert to Z sores first
CTzTrials = CT(:,:,zTrials);

CTZ = nanzscore(reshape(CTzTrials,...
    size(CTzTrials,1),size(CTzTrials,2)*size(CTzTrials,3))')';
CTZ = reshape(CTZ, size(CTzTrials,1), size(CTzTrials,2), size(CTzTrials,3));

% relabel trials
zObjectTrials = repmat([1 0 1 0 1 0], 1, 8);
zSocialTrials = repmat([0 1 0 1 0 1], 1, 8);
zStimTrials = repmat([1 1 1 1 1 1 0 0 0 0 0 0], 1, 4);

% K MEANS CLUSTERING OF Z SCORES %%

CZS = CTZ(:,:,find(zStimTrials));
CZE = CTZ(:,:,find(abs(zStimTrials)-1));

%
CZOS = CZS(:,:,1:2:end);
CZSS = CZS(:,:,2:2:end);
CZOE = CZE(:,:,1:2:end);
CZSE = CZE(:,:,2:2:end);

X = {CZOS CZOE CZSS CZSE};
for i = 1:numel(X)
    X{i} = nanmean(X{i},3);
end


%%
close all
kClusters = 5; %define the number of k means pairs
titles = {'Object with Opto' 'Object No Opto' 'Social with Opto' 'Social No Opto'};
        
%xaxis
        %Make X Axis Correctly
        xlabel('Time (s)')
        x = (1:size(X{1},2))*(1/30.98);
        x = x(1:30.98*5:end);
        xticks = 1:size(X{1},2);
        xticks = xticks(1:31*5:end);


for ii = [ 1 2 3 4] % for each trial type
    
    %perform k means clustering
    tic
    [IDX, C, SUMD, D] = kmeans(X{ii}, kClusters);
    toc

    for i = 1:numel(X)
        Y = [];
        for j = 1:kClusters
            Y = [Y; X{i}(find(IDX==j),:)];
            YY{ii,j,i} =  X{i}(find(IDX==j),:);

            
        end
        XClustered{i} = Y;   
    end
    
    %figure
    h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 5 2];
    for i = 1:numel(X)
        subplot(1,numel(X),i)
        imagesc(XClustered{i},[0 2])
        colormap(flipud(gray))
        
        %title
        title(titles{i})   
        
        %yaxis
        if i == 1
            ylabel('Cortical Neuron Number')
        else       
            set(gca,'YTick',[])
        end
        
       
        set(gca,'XTick', xticks, 'XTickLabel',round(x), 'FontSize',3);
 
       
        PlotVerticalLines([...
            nFramesBaseline ...
            nFramesBaseline + nFramesAdvance ...
            nFramesBaseline + nFramesAdvance + nFramesStationary...
            nFramesBaseline + nFramesAdvance + nFramesStationary + nFramesRetreat], 0 , 2000)
        
        PlotVerticalLines([...
            nFramesBaseline ...
            nFramesBaseline + nFramesAdvance ...
            nFramesBaseline + nFramesAdvance + nFramesStationary ...
            nFramesBaseline + nFramesAdvance + nFramesStationary + nFramesRetreat],  12000, 13800)
       
    end
    saveas(h101,fullfile(fig_dir, titles{ii}),'pdf')
    
    
end

%save('~/2presults/8armStim/allNeuronsTrialZscores.mat', '-v7.3')

%
c = {'r' 'm' 'b' 'c'}
for i = 1:size(YY,1)
    h102 = figure; h102.Units = 'inches'; h102.Position = [1 1 2 5];
    
    for j = 1:size(YY,2)
        hold on;
        subplot(5,1,j)
        
        for k = 1:size(YY,3)
            hold on
            shadedErrorBar([],mean(YY{i,j,k}), sem(YY{i,j,k}),c(k))   
        end
        
        %ylabel
        ylabel 'ZScore'
        
        %Make X Axis Correctly
        xlabel('Time (s)')
        x = (1:size(YY{1,1,1},2))*(1/30.98);
        x = x(1:31*5:end);
        xticks = 1:size(YY{1,1,1},2);
        xticks = xticks(1:31*5:end);
        set(gca,'XTick', xticks, 'XTickLabel',round(x), 'FontSize',5);
        axis([0 1082 -inf inf])

    end
    %saveas(h102,fullfile(fig_dir, [titles{i} 'KMeansTraces']),'pdf')

end

        
  %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PCA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%XXX = [X{1} X{2} X{3} X{4}];
nFramesTrial = round(nFramesAdvance+nFramesStationary+nFramesRetreat);
XXX = [X{1}(:,nFramesBaseline:nFramesBaseline+nFramesTrial) ...
    X{2}(:,nFramesBaseline:nFramesBaseline+nFramesTrial)...
    X{3}(:,nFramesBaseline:nFramesBaseline+nFramesTrial) ...
    X{4}(:,nFramesBaseline:nFramesBaseline+nFramesTrial)];
if ~exist('COEFF')
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(XXX);
end
figure; hold on
%t = [1 1082 1082*2 1082*3 1082*4]
t = [1 nFramesTrial+2 nFramesTrial*2+3 nFramesTrial*3+4 nFramesTrial*4+4]
for i = 1:4
    h = scatter3(...
            COEFF(t(i) : t(i+1),1)-COEFF(t(i),1), ...
            COEFF(t(i) : t(i+1),2) -COEFF(t(i),2),...
            COEFF(t(i) : t(i+1),3) - COEFF(t(i),3),...
            c{i},'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha', 1)
        
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')

end
%%

    for i = 1:5
        figure
        for j = 1:4
    hold on
    plot(COEFF(t(j) : t(j+1),i),c{j})

end
end



%%



%get one value for each trial from each neuron
CTZA = [];
for i = 1:64
    CTZA(:,i) = nanmean(CT(:, stimStartFrame:stimEndFrame,i),2);
end


%%
COZS = (CTZA .* stim) .* objectTrials;
COZE = (CTZA .*abs(stim -1)) .* objectTrials;

COZS(:, all(~any(COZS),1)) = [];
COZE(:, all(~any(COZE),1)) = [];


%calcuate a p value for stim vs nostim object trials
p = ones(size(COZS,1),1);
for i = 1:size(COZS,1)
    p(i) = ranksum(COZS(i,:), COZE(i,:));
end












%% INDIVIDUAL NEURON ANLAYSIS --STIM NO STIM -- USING SPIKING DATA 
%how does the a neuron's activity change with stimulation

% % %for OBJECT trials
SOS = (STS .* stim) .* objectTrials;
SOE = (STS .*abs(stim -1)) .* objectTrials;

SOS(:, all(~any(SOS),1)) = [];
SOE(:, all(~any(SOE),1)) = [];

%calcuate a p value for stim vs nostim object trials
p = ones(size(SIS,1),1);
for i = 1:size(SOS,1)
    p(i) = ranksum(SOS(i,:), SOE(i,:));
end

%find the cells that are significantly different
I_soc = find(p < 0.01)


% % % for SOCIAL trials
SSS = (STS .* stim) .* socialTrials;
SSE = (STS .*abs(stim -1)) .* socialTrials;

SSS(:, all(~any(SSS),1)) = [];
SSE(:, all(~any(SSE),1)) = [];

%calcuate a p value for stim vs nostim object trials
p = ones(size(SSS,1),1);
for i = 1:size(SOS,1)
    p(i) = ranksum(SSS(i,:), SSE(i,:));
end

%find the cells that are significantly different
I_soc = find(p < 0.01)


%% Let's plot all of these cells that are different.
plotObjectCells = 1
if plotObjectCells
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 8 2]
hold on
for i = 1:numel(I_obj)
    scatter((i)*ones(1,12), SOS(I_obj(i),:),'r','filled')
    scatter((i+0.25)*ones(1,12),  SOE(I_obj(i),:),'k','filled')
end

xlabel 'Neuron Numbr'
ylabel 'Average Fl'
end


plotSocialCells = 1
if plotSocialCells
h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 8 2]
hold on
for i = 1:numel(I_soc)
    scatter((i)*ones(1,12), SSS(I_obj(i),:),'r','filled')
    scatter((i+0.25)*ones(1,12),  SSE(I_obj(i),:),'k','filled')
end

xlabel 'Neuron Numbr'
ylabel 'Average Fl'
end

%plot a venn diagram of significantly modulated neurons
plotVennDiagram = 1
if plotVennDiagram
    I_intersect = intersect(I_obj, I_soc);
    h3 = figure; h3.Units = 'inches'; h3.Position = [5 5 2 1]
    venn([length(I_obj) length(I_soc)], length(I_intersect))
    saveas(h3,fullfile(fig_dir, 'Venn'),'pdf')
end

%% INDIVIDUAL NEURON ANALYSIS --OBJ VS SOC NO STIM -- 
COE = (CTA .*abs(stim -1)) .* objectTrials;
CSE = (CTA .*abs(stim -1)) .* socialTrials;

%calcuate a p value for object vs social trials
p = ones(size(COE,1),1);
for i = 1:size(COE,1)
    p(i) = ranksum(COE(i,:), CSE(i,:));
end

%find the cells that are significantly different
I = find(p < 1)

plotVennDiagram = 1
if plotVennDiagram
    I_intersect = intersect(I_obj, I_soc);
    h3 = figure; h3.Units = 'inches'; h3.Position = [5 5 2 1]
    venn([length(I_obj) length(I_soc)], length(I_intersect))
    saveas(h3,fullfile(fig_dir, 'Venn'),'pdf')
end

%scatter figure
scatter(nanmean(COE'), nanmean(CSE'),'filled')


%% POPULATION LEVEL ANALYSIS








%% CLASSIFIERS AND MODELS


%KMeans








