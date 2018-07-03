 
load('~/2presults/8arm2/allNeuronsTrial.mat')

%Important for this experiment
nFramesStationary = 212;
nFramesBeforeStationary = nFramesBeforeStim + (nFramesStim - nFramesStationary);

%INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
cTO = CT(:,:,find(objectTrials));
cTS = CT(:,:,find(socialTrials));
cTB = CT(:,:,find(blankTrials));
cTT = CT(:,:,find(toyTrials));

d = {cTO cTS cTT cTB};
labels = {'object' 'social' 'toy' 'blank'}


%Find Peaks and Locations of Calcium Traces

PKS = {}; LOCS = {}; 

% count the number of calcium spikes
for i = 1:numel(d)
    for j = 1:size(CT,1)
 
       % extract calcium peak times
       for k = 1:size(d{i},3)
           [PKS{i,j,k}(:)] = findpeaks(d{i}(j,:,k));
           [~, LOCS{i,j,k}(:)] = findpeaks(d{i}(j,:,k)*-1);
       end
        
        %extract spikes peak times  
%         for k = 1:size(d{i},3)
%             [PKS{i,j,k}(:)] = sum(d{i}(j,:,k));
%             [~, LOCS{i,j,k}(:)] = findpeaks(d{i}(j,:,k));
%         end
        
        
        
    end
end

%count spikes
C_spikes = []; 
for i = 1:size(CT,1) 
        tP = [];
        
        %object
        x = 0;
        for k = 1:size(cTO,3)
           tP(k,1) = numel(LOCS{1,i,k});
            x = x + numel(LOCS{1,i,k})
        end
        C_spikes(i,1) = x;
        
        
        %social
        y = 0;
        for k = 1:size(cTO,3)
            tP(k,2) = numel(LOCS{2,i,k});
            y = y + numel(LOCS{2,i,k})
        end
        
        
        C_spikes(i,2) = y;
        
        %obj vs social
        C_spikes(i,3) = ranksum(tP(:,1), tP(:,2));

        
        %novelty
        C_spikes(i,4) = ranksum([tP(1:5,1) ; tP(1:5,2)],[tP(19:24,1);tP(19:24,2)]);

end

%%
%[p, index] = sort(C_spikes(:,3))

%Lets plot some examples
plotTraces = 1
if plotTraces 
for i = 270
    
    h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 1 3]; 
    
    subplot(3,1,1);
    
    
    
    plot(mean(cTO(i,:,:),3),'r') ; hold on
    plot(mean(cTS(i,:,:),3),'b')
    %plot(mean(cTB(index(i),:,:),3),'k')
    %plot(mean(cTT(index(i),:,:),3),'g')
    
    ylabel('Fluorescence (trial avg.)')
     %Make X Axis Correctly
        xlabel('Time (s)')
        x = (1:size(CT,2))*(1/31);
        x = round(x(1:31*5:end));
        xticks = 1:size(CT,2);
        xticks = xticks(1:31*5:end);
        title('Average Calcium Traces')
        
          PlotVerticalLines([nFramesBeforeStationary nFramesBeforeStimEnd])
         set(gca,'XTick',xticks,'XTickLabel',x,'FontSize',5)
            axis([-0 713 0 7])

         %object
        subplot(3,1,2)
       for k = 1:size(LOCS,3)
               hold on
               plot([LOCS{1,i,k}; LOCS{1,i,k}], [k*ones(length(LOCS{1,i,k}),1)-.4, k*ones(length(LOCS{1,i,k}),1)+0.4]','k','LineWidth',.3) 
               
       end
       set(gca,'XTick',xticks,'XTickLabel',x)
       xlabel('Times (s)')
       ylabel('Trial Number')
         PlotVerticalLines([nFramesBeforeStationary nFramesBeforeStimEnd])
         set(gca,'XTick',xticks,'XTickLabel',x,'FontSize',5)
        axis([-0 713 0 25])
        title('Object Trials')

       %social
        subplot(3,1,3)
        title('Social Trials')
       for k = 1:size(LOCS,3)
               hold on
               plot([LOCS{2,i,k}; LOCS{2,i,k}], [k*ones(length(LOCS{2,i,k}),1)-.4, k*ones(length(LOCS{2,i,k}),1)+0.4]','k','LineWidth',.3) 
               
       end
       set(gca,'XTick',xticks,'XTickLabel',x)
       xlabel('Times (s)')
       ylabel('Trial Number')
         PlotVerticalLines([nFramesBeforeStationary nFramesBeforeStimEnd])
         set(gca,'XTick',xticks,'XTickLabel',x,'FontSize',5)
        axis([0 713 0 25])

                saveas(h1, ['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm2/Figures/SingleNeuron/'  num2str(i)],'pdf')
        
end
end





%% %% %% %% %% %
%Further analysis

%neuron 270 
num = 8
i(num)
p1 = []; p2 = [];
for k = 1:size(PKS,3)
    p1(k) = numel(PKS{1,i(num),k});
    p2(k) = numel(PKS{2,i(num),k});
end

h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 1 1];

hold on
histogram(p1,'FaceColor','r','EdgeAlpha',0.1,'EdgeColor','r')
histogram(p2,'FaceColor','b','EdgeAlpha',0.1,'EdgeColor','b')
ylabel('Num of Trials')
xlabel('Num Spikes / Trial')
set(gca,'FontSize',5)


saveas(h2, ['~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm2/Figures/SingleNeuron/Histogram_' num2str(i(num))],'pdf')
