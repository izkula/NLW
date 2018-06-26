

clear; clc; close all
tic
load('~/2presults/8armStim/allNeuronsTrial.mat')
toc




%%
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
labels = {'object' 'social' 'obj stim' 'soc stim'}



%%%%%%%%%%%%%%%%%%%%%%%% KMEANS CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%put the organized z scores into one cell array
X = {CZO CZOS CZS CZSS};
for i = 1:numel(X)
    X{i} = nanmean(X{i},3);
end


%% Save Load Point
clear CTZ CTzTrials ST  AA ans bpodDataPath bpodImagePath
%save('~/2presults/8armStim/Z.mat','-v7.3')
load('~/2presults/8armStim/Z.mat')

%%  K MEANS CLUSTERING
close all
kClusters = 4; %define the number of k means pairs
titles = {'O' 'OS' 'S' 'SS'  };
c = flipud(linspecer(kClusters, 'distinguishable'))

%xaxis
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:size(X{1},2))*(1/30.98);
x = round(x(1:30.98*5:end));
xticks = 1:size(X{1},2);
xticks = xticks(1:31*4:end);

C = {}; C_trace_avg = {}
for ii = 1:4 % for each trial type
    %perform k means clustering
    tic
    [IDX, C{ii}, SUMD, D] = kmeans(X{ii}, kClusters);
    toc
      
    doTracePlot = 1
    if doTracePlot
    h102 = figure; h102.Units = 'inches'; h102.Position = [1 1 3.5 1];

    for i = 1:numel(X)
        Y = [];
        
        for j = 1:kClusters
            Y = [Y; X{i}(find(IDX==j),:)];  
            pause(0.1)
            %generate a plot for the mean of each k groups
            subplot(1,numel(X),i); hold on
            plot(mean(X{i}(find(IDX==j),:)),'Color',c(j,:))  
            if ii == i
                C_trace_avg{ii,j} =  mean(X{i}(find(IDX==j),:))   
            end 
        end
        
        %some plotting features and saving
        ylabel('Z Score')
        xlabel('time (s)')
        set(gca,'XTick',xticks,'XTickLabels',x,'FontSize',3)
        axis([-inf inf -1 1])
        
        saveas(h102, fullfile(fig_dir, 'Kmeans', [titles{ii} '_traces']),'pdf')
        %save out to a cell array
        XClustered{i} = Y;   
    end
    
    end
    %figure
    
    
    doBigPlot = 1
    if doBigPlot 
    h101 = figure; h101.Units = 'inches'; h101.Position = [1 1 3.5 1.5];
   
    for i = 1:numel(X)
        subplot(1,numel(X),i)
        imagesc(XClustered{i},[0 1.5])
        %colormap(flipud(gray))
        
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
             nFramesBeforeStim, ...
             nFramesBeforeStim + nFramesStim])
       
    end
    saveas(h101,fullfile(fig_dir, 'Kmeans', titles{ii}),'pdf')
    close  
    
    end

end

%% Let's do some statistics on the trial averaged data
   %for each neuron, let's find the correlation of the trial averaged trace
   %for the four conditions
   R = [];
   for i = 1:size(X{1},1)
       R(:,:,i) = corr([X{1}(i,:)' X{2}(i,:)' X{3}(i,:)' X{4}(i,:)']);     
   end
   
Ravg = mean(R,3);

%figure
h103 = figure; h103.Units = 'inches'; h103.Position = [1 1 1.75 1];
imagesc(Ravg, [-0 1]); colorbar; 


set(gca,'XTick',[1:4],'XTickLabel',{'Object' 'Social' 'Object+Stim' 'Social+Stim'},...
    'YTick',[1:4],'YTickLabel',{'Object' 'Social' 'Object+Stim' 'Social+Stim'},'FontSize',3,'XTickLabelRotation',45)

saveas(h103, fullfile(fig_dir, 'Kmeans','CorrelationsSummary'),'pdf')

 
%% Now let's calculate the correlation and the k means  for each mouse
cellCounterEachMouse = [1; cumsum(sum(cellCounter,2))]

XbyMouse = {}
for i = 1:length(cellCounter)
    for j = 1:numel(X)
    XbyMouse{i,j} = X{j}(cellCounterEachMouse(i):cellCounterEachMouse(i+1),:)
    
    end   
end
%%
    




R_eachMouse = [];
for i = 1:length(cellCounterEachMouse)-1
    x = R(:,:,cellCounterEachMouse(i):cellCounterEachMouse(i+1));
    R_eachMouse(:,:,i) = mean(x,3)
    
end

R_eachMouse = reshape(R_eachMouse,[],8)';
R_eachMouse = R_eachMouse(:,[2 3 4 7 8  12  ])
R_eachMouse_avg = mean(R_eachMouse);
R_eachMouse_sem = sem(R_eachMouse);

h104 = figure; h104.Units = 'inches'; h104.Position = [1 1 3 1];
hold on
%(R_eachMouse_avg,R_eachMouse_sem)
bar(R_eachMouse_avg);
for i = 1:6
    scatter(i*ones(1,8),R_eachMouse(:,i),1.5,'MarkerFaceColor','k','MarkerEdgeColor','k')
end

ylabel('Pearsons Corr Each Mouse')
set(gca,'XTick',1:9,'FontSize',3,'XTickLabelRotation',45,'XTickLabels',{...
    'Object-Social', 'Object-ObjectStim', 'Object-SocialStim',...
    'Social-ObjectStim','Social-SocialStim',...
    'ObjectStim-SocialStim'})
axis([-inf inf -0.5 0.5])
saveas(h104, fullfile(fig_dir,'Kmeans','CorrelationsEachMouse'),'pdf')
    
%% Now Let's calculate teh K means clusters for each mouse independently

kmeans_single_mouse(XbyMouse, cellCounter,5)


  %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PCA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
close all
XXX = [X{1} X{2} X{3} X{4}];
%XXX = [X{1}(:,nFramesAdvance-31:nFramesAdvance+nFramesStationary+31) ...
%    X{2}(:,nFramesAdvance-31:nFramesAdvance+nFramesStationary+31)...
%    X{3}(:,nFramesAdvance-31:nFramesAdvance+nFramesStationary+31) ...
%    X{4}(:,nFramesAdvance-31:nFramesAdvance+nFramesStationary+31)];

%if ~exist('COEFF')
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(XXX);
%end
x = size(COEFF,1)/4
C2 = mat2cell(COEFF,[x x x x] , size(COEFF,2))

for i = 1:4
    for j = 1:size(C2,2)
        C2{i}(:,j) = smooth(C2{i}(:,j),62);
    end
end



%% Make a Plot of the PCS
figure; hold on
c = {'r' 'm' 'b' 'c'}
for i = 1:4
    
    %3d plot
    h = plot3(...
        C2{i}(:,1) - C2{i}(1,1) , ...
        C2{i}(:,2) - C2{i}(1,2),...
        C2{i}(:,3) - C2{i}(1,3), c{i},'LineWidth',1) 
        
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    plot3(C2{i}(1,1),C2{i}(1,2), C2{i}(1,3), 'k*')
   % plot3(C2{i}(nFramesBeforeStim,1),C2{i}(nFrames,2), C2{i}(nFramesAdvance,3), 'ks')
   % plot3(C2{i}(nFramesAdvance+nFramesStationary,1),...
   %    C2{i}(nFramesAdvance+nFramesStationary,2), C2{i}(nFramesAdvance+nFramesStationary,3), 'ko')
   % plot3(C2{i}(end-nFramesPostTrial,1),C2{i}(end-nFramesPostTrial,2), C2{i}(end-nFramesPostTrial,3), 'kd')
   plot3(C2{i}(end,1),C2{i}(end,2), C2{i}(end,3), 'k+')
    
end
saveas(h,[fig_dir 'PC3d'],'fig')


for i = 1:10
    h1 = figure; hold on
    %each PC
    title(num2str(i))
    for j = 1:numel(C2)
     plot(C2{j}(:,i),c{j})   
    end
    
    
end


%% Calculate Distance

nPCs = 1;
D_stim = nan(1,928); D_nostim = nan(1,928);
for i = 1:size(C2{1},1)
    D_stim(i) = norm(C2{1}(i,1:nPCs) - C2{2}(i,1:nPCs))  ;
    D_nostim(i) = norm(C2{3}(i,1:nPCs) - C2{4}(i,1:nPCs));
end
%%
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 3 2];
plot(D_stim,'b')
hold on
plot(D_nostim,'k')

PlotVerticalLines([nFramesAdvance nFramesAdvance+nFramesStationary nFramesAdvance+nFramesStationary+nFramesRetreat])
%ylabel
ylabel 'Distance (Object - Social , a.u.)'
        
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:nFramesTrial)*(1/30.98);
x = x(1:31*5:end);
xticks = 1:nFramesTrial;
xticks = xticks(1:31*5:end);
set(gca,'XTick', xticks, 'XTickLabel',round(x), 'FontSize',5);
axis([0 928 -inf inf])
saveas(h1,[fig_dir 'PC_Distance'],'pdf')

%%

for i = 1:5
    h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 1.5 2];
  for j = 1:4
    hold on
    plot(C2{j}(:,i), c{j})

end
end
PlotVerticalLines([nFramesAdvance nFramesAdvance+nFramesStationary nFramesAdvance+nFramesStationary+nFramesRetreat])
%ylabel
ylabel 'PC 3'
        
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:nFramesTrial)*(1/30.98);
x = x(1:31*5:end);
xticks = 1:nFramesTrial;
xticks = xticks(1:31*5:end);
set(gca,'XTick', xticks, 'XTickLabel',round(x), 'FontSize',5);
axis([0 928 -inf inf])
saveas(h2,[fig_dir 'PC3'],'pdf')







%% Cosine Distnace


nPCs = 4;
D_stim = nan(1,nFramesTrial); D_nostim = nan(1,nFramesTrial);
for i = 1:size(C2{1},1)
    D_stim(i) = pdist2(C2{1}(i,1:nPCs),C2{3}(i,1:nPCs),'cosine')  ;
    D_nostim(i) = pdist2(C2{2}(i,1:nPCs), C2{4}(i,1:nPCs),'cosine');
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








