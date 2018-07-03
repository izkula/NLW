%clear; clc; close all
tic
load('~/2presults/8armStim/allNeuronsTrial.mat')
toc
%

%%Calculate Z Score Across Entire Data Set
CT_reshaped = reshape(CT,...l
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
%%
r = size(CT,1); col = 50;

d_m = [mean(d{1},3) nan(r,col)  mean(d{2},3) nan(r,col) mean(d{3},3) nan(r,col) mean(d{4},3)];

kClusters = 4
c = flipud(linspecer(kClusters, 'distinguishable'))

d_all = {}
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 11 8];
for i = 1:4
    [IDX, C, SUMD, D] = kmeans(mean(d{i},3), kClusters);

    d_m_sorted = [];
     for j = 1:kClusters
         d_m_sorted = [d_m_sorted;  d_m(find(IDX==j),:)];
         d_all{i,j} = d_m(find(IDX==j),:);

     end
     
     
     doPlotAll = 1
     if doPlotAll 
         %plot
         subplot(1,4,i)
         imagesc(d_m_sorted,[0 1.5])

         drawFrames = [nFramesBeforeStim nFramesBeforeStationary nFramesBeforeStim + nFramesStim]
         
         PlotVerticalLines([(0*(col+nFramesTrial)+ drawFrames)  ])
         PlotVerticalLines([(1*(col+nFramesTrial)+ drawFrames)  ])
         PlotVerticalLines([(2*(col+nFramesTrial)+ drawFrames)  ])
         PlotVerticalLines([(3*(col+nFramesTrial)+ drawFrames)  ])

         ylabel('Neuron ID #')
         xlabel('Time (s)')


            %xaxis
            %Make X Axis Correctly
            xlabel('Time (s)')
            x = (1:size(CT,2))*(1/30.98);
            x = round(x(1:30.98*4:end));
            xticks = 1:size(CT,2);
            xticks = xticks(1:31*4:end);
   
            
            
            
            set(gca,'XTick',xticks,'XTickLabels',x,'FontSize',5)
            if i == 4
               fig_dir = '~/Dropbox/2p_Claustrum_Shared/2p/Results/8armStim/Figures/'
               saveas(h1,[fig_dir '/KMeansAll'], 'pdf')
            end
     end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lets find the footprint of each neuron, grouped by kmeans cluster%

cellIdx = {}; counter = 0; A_K = {}
for i = 1:size(cellCounter,1)
    for j = 1:size(cellCounter,2)   
        x = [];
        
        cellIdx{i,j} = IDX(counter+1:counter+cellCounter(i,j))
        counter = counter + cellCounter(i,j)
        
        
        
        
        
        
        
        
        
        
        
        
        %x = full(AA{i,j});
        
        %could make 0 and 1 here
        %x(x>0)=1;
        
        %y = reshape(x,254,300,size(x,2));
        
        
        %calculate the center 
        
        
        
        %for k = 1:kClusters
            
            
            
            
            
            
            
        %end
        %
        
        
        
        
    end
end
    
    

    







