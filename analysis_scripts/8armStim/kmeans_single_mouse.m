function [outputArg1,outputArg2] = kmeans_single_mouse(X,cellCounter, kClusters)
%colors
c = flipud(linspecer(kClusters,'distinguishable'))

%xaxis
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:size(X{1},2))*(1/30.98);
x = round(x(1:30.98*5:end));
xticks = 1:size(X{1},2);
xticks = xticks(1:31*5:end);


for ii = 1:length(cellCounter)
    for jj = 1:size(X,2) 
        ii
        jj
        [IDX, C, SUMD, D] = kmeans(X{ii,jj},kClusters); % k means clustering

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        doTracePlot = 0
        if doTracePlot
            h102 = figure; h102.Units = 'inches'; h102.Position = [1 1 3.5 1]; %make a figure for each mouse x each trial type

            for i = 1:size(X,2)
                Y = [];

                for j = 1:kClusters
                    Y = [Y; X{ii,jj}(find(IDX==j),:)];  
                    pause(0.1)
                    %generate a plot for the mean of each k groups
                    subplot(1,size(X,2),i); hold on
                    plot(mean(X{ii,i}(find(IDX==j),:)),'Color',c(j,:))  
                    %if ii == i
                    %    C_trace_avg{ii,j} =  mean(X{i}(find(IDX==j),:))   
                    %end 
                end

                %some plotting features and saving
                ylabel('Z Score')
                xlabel('time (s)')
                set(gca,'XTick',xticks,'XTickLabels',x,'FontSize',3)
                axis([-inf inf -1 1])
                
                XClustered{i,j} = Y;   

            end
                %save out to a cell array
                %saveas(h102, fullfile(fig_dir, 'Kmeans', [titles{ii} '_traces']),'pdf')

        end
    end
end


end
