
clear; clc; close all
tic
load('~/2presults/8armStim/allNeuronsTrial.mat')
toc



%INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
cTO = CT(:,:,find(objectTrials.*abs(stimTrials-1)));
cTS = CT(:,:,find(socialTrials.*abs(stimTrials-1)));
cTOS = CT(:,:,find(objectTrials.*stimTrials));
cTSS = CT(:,:,find(socialTrials.*stimTrials));

d = {cTO cTOS cTS cTSS};
labels = {'object' 'social' 'obj stim' 'soc stim'}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PCA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
framesUsedForPCA = 1 : nFramesTrial;
plotOffset = 0 ;
%calculate means
d_m = [mean(d{1}(:,framesUsedForPCA,:),3)...
    mean(d{2}(:,framesUsedForPCA,:),3) mean(d{3}(:,framesUsedForPCA,:),3)...
    mean(d{4}(:,framesUsedForPCA,:),3)];


d_m_z = zscore(d_m')';

%run PCA
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(d_m_z);

x = size(d_m,2)/4
C2 = mat2cell(COEFF,[x x x x] , size(COEFF,2))


%smooth and center
for i = 1:4
    for j = 1:size(C2{1},2)
        C2{i}(:,j) = smooth(C2{i}(:,j) - C2{i}(1,j)  ,  31) ;
    end
end


% Make a Plot of the PCS
c = {'r' 'm' 'b' 'c'}

fig_dir = '~/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
plotPCA = 0
if plotPCA
    
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 2.5 2];
hold on
for i = 1:4 
    %3d plot
    plot3(...
        C2{i}(:,1) , ...
        C2{i}(:,2),...
        C2{i}(:,3), c{i},'LineWidth',1) 
        
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    plot3(C2{i}(1,1),C2{i}(1,2), C2{i}(1,3), 'k*')
    plot3(C2{i}(end,1),C2{i}(end,2), C2{i}(end,3), 'k+')
    plot3(C2{i}(end-nFramesAfterStimEnd-plotOffset,1),C2{i}(end-nFramesAfterStimEnd-plotOffset,2), ...
        C2{i}(end-nFramesAfterStimEnd-plotOffset,3), 'kd')

if i == 2 || i == 4
    plot3(C2{i}(nFramesBeforeStim-plotOffset,1),C2{i}(nFramesBeforeStim-plotOffset,2), C2{i}(nFramesBeforeStim-plotOffset,3), 'ks')

end
end
set(gca,'FontSize', 5)
saveas(h1,[fig_dir 'PC3d'],'pdf')
end
%%
close all
for i = 1:4
    h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 2.5 1.5]; 
    hold on
    %each PC
    title(num2str(i))

    for j = 1:numel(C2)
     plot(C2{j}(:,i),c{j})  

     ylabel(['PC ' num2str(i) ' a.u.'])
     
    end
        nFramesStationaryEnd = size(CT,2) - nFramesAfterStimEnd;

    plot([nFramesBeforeStim nFramesBeforeStim], [-0.005 0.005],'k')
    %plot([nFramesBeforeStationary nFramesBeforeStationary], [-0.005 0.005],'k')
    plot([nFramesStationaryEnd nFramesStationaryEnd], [-0.005 0.005],'k')
    
         %Make X Axis Correctly
    xlabel('Time (s)')
    x = (1:size(CT,2))*(1/31)
    x = round(x(1:31*4:end))
    xticks = 1:size(CT,2)
    xticks = xticks(1:31*4:end)
    
    axis([0 713 -inf inf])
    set(gca,'XTick',xticks,'XTickLabel',x)

    %lines
    saveas(h1, [fig_dir num2str(i)], 'pdf')

end
%%
%Supplemental Subplotting
    h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 5 4]; 

for i = 1:10
    subplot(5,2,i)
    %each PC
    title(num2str(i))

    for j = 1:numel(C2)
        hold on
     plot(C2{j}(:,i),c{j})  

     ylabel(['PC ' num2str(i) ' a.u.'])
     
     %Make X Axis Correctly
    xlabel('Time (s)')
    x = (1:size(CT,2))*(1/31);
    x = round(x(1:31*4:end));
    xticks = 1:size(CT,2);
    xticks = xticks(1:31*4:end);
    
    axis([0 713 -inf inf])

    end
        nFramesStationaryEnd = size(CT,2) - nFramesAfterStimEnd;

    plot([nFramesBeforeStim nFramesBeforeStim], [-0.005 0.005],'k')
   % plot([nFramesBeforeStationary nFramesBeforeStationary], [-0.005 0.005],'k')
    plot([nFramesStationaryEnd nFramesStationaryEnd], [-0.005 0.005],'k')
    
    set(gca,'XTick',xticks,'XTickLabel',x,'FontSize', 3)

    %lines

end
saveas(h1, [fig_dir 'PCs20'], 'pdf')


%%
%Supplemental Explained Variance

h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 2 2];
plot(cumsum(EXPLAINED),'ok')
axis([0 10 0 100])
ylabel('Cum. % Explained')
xlabel(' PC #')
set(gca, 'FontSize',5)
saveas(h2, [fig_dir 'ExplainedCA'],'pdf')


%% PCA each mouse AND EUCLIDEAN DISTANCE %%%

cellNumbers = [1 cumsum(sum(cellCounter')) length(cellCounter)]

dist = {};
for ii = 1:size(cellCounter,1)
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(d_m_z(cellNumbers(ii):cellNumbers(ii+1),:));


x = size(d_m,2)/4
C2 = mat2cell(COEFF,[x x x x] , size(COEFF,2))

    
%smooth and center
for i = 1:4
    for j = 1:size(C2{1},2)
        C2{i}(:,j) = smooth(C2{i}(:,j) - C2{i}(1,j),  31) ;
    end
end


% Make a Plot of the PCS
c = {'r' 'm' 'b' 'c'}

fig_dir = '~/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
plotPCA = 0
if plotPCA
    
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 3 2];
subplot(1,6,1:4)
hold on
for i = 1:4 
    %3d plot
    plot3(...
        C2{i}(:,1) , ...
        C2{i}(:,2),...
        C2{i}(:,3), c{i},'LineWidth',1) 
        
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    plot3(C2{i}(1,1),C2{i}(1,2), C2{i}(1,3), 'k*')
    plot3(C2{i}(end,1),C2{i}(end,2), C2{i}(end,3), 'k+')
    plot3(C2{i}(end-nFramesAfterStimEnd-plotOffset,1),C2{i}(end-nFramesAfterStimEnd-plotOffset,2), ...
        C2{i}(end-nFramesAfterStimEnd-plotOffset,3), 'kd')

if i == 2 || i == 4
    plot3(C2{i}(nFramesBeforeStim-plotOffset,1),C2{i}(nFramesBeforeStim-plotOffset,2), C2{i}(nFramesBeforeStim-plotOffset,3), 'ks')

end
set(gca,'FontSize',5)
end
subplot(1,6,5:6)
plot(cumsum(EXPLAINED),'o','MarkerSize',5)
axis([0 10 0 100])

set(gca,'FontSize', 5)
saveas(h1,[fig_dir 'PC3d_Mouse' num2str(ii)],'pdf')


end

calculatePCDistance = 1;
if calculatePCDistance
    nPCsToUse = 4;
    
    d_proj = COEFF(:,1:nPCsToUse) * SCORE(:,1:nPCsToUse)';
  
    x = size(d_proj,1)/4
    C2 = mat2cell(d_proj,[x x x x] , size(d_proj,2));
    
%calculate distancesfor each mouse
    IDXs = [ 1 3; 1 2; 3 4; 2 4; ]
    for j = 1:size(IDXs,1)
        for i = 1:x
            
            u = C2{IDXs(j,1)}(i,:);
            v = C2{IDXs(j,2)}(i,:);
            
            %dist{j}(i,ii) = pdist(u,v)
            %dist{j}(i,ii) =  dist{j}(i,ii) ;
             
            dist{j}(i,ii) = dot(u,v)/(norm(u)*norm(v));
             
             
        end     
    end 
end


end

%%    EUCLIDEAN DISTANCE
  % plotting
 h4 = figure; h4.Units = 'inches'; h4.Position = [1 1 4.5 1.25]

titles = {'Object / Social' 'Object / Object with Stim' 'Social / Social with Stim' 'Object with Stim / Social with Stim'}
for i = 1:numel(dist)
    subplot(1,4,i)

    hold on
    temp = dist{i} - dist{i}(1,:)


    ShadedError(1:713,mean(temp'),sem(temp'))


    plot(nFramesBeforeStim,mean(temp(nFramesBeforeStim,:)),'rs')
    plot(nFramesBeforeStim+nFramesStim,mean(temp(nFramesBeforeStim+nFramesStim,:)),'rd')
    
     %Make X Axis Correctly
    xlabel('Time (s)')
    x = (1:size(CT,2))*(1/31);
    x = round(x(1:31*4:end));
    xticks = 1:size(CT,2);
    xticks = xticks(1:31*4:end);
    set(gca,'XTick',xticks,'XTickLabels',x,'FontSize',5)
    axis([-inf inf -inf inf])
    if i == 1
    ylabel('Euclidean Distance')
    end
saveas(h4, [fig_dir 'EuDist' num2str(i)],'pdf')
end 






