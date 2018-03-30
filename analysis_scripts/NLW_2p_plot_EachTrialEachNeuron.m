clear; close all; clc
init_samz
%%The script imports trialtraces and then plots all neurons%%

%Logic
doIndividualStimulusFigures = 1

f = dir('~/2presults/8arm/')

%aggregates trace data from all mice and trials and creates one trial
%averaged matrix

trialsOrdered = [1 2 8 3 5 4 6 7]
titles = {'Obj 1' 'Soc 1' 'Empty' 'Obj 2'  'Soc 2' 'ToySoc' 'Obj3' 'Soc 3'} 
neuronNumber = 1
for i = 3:numel(f)
    dirName = fullfile(f(i).folder,f(i).name);
    
    %load neuron trial trace data
    load(fullfile(dirName, 'trialtraces.mat'))
    
    for j = 1:size(traces{1},1) %for each neuron
        try
       h = figure; h.Units = 'inches'; h.Position = [1 1 8 8];
       set(gcf, 'Visible','off')
       
        for k = 1:numel(trialsOrdered) %for each stimulus ordered by type
           
            %for l = 1:size(traces{1},3) %for each trial
             %   hold on
              %  d = traces{k}(j,:,l);
                                
                %subplot(3,3,k)
               
                %t = linspace(-15,length(d)/frameSampleRate-15,length(d));

                %plot(t,d,'k')        

            %end
            
            %shadedErrorBar(t,mean(traces{k}(j,:,:),3),...
                %sem(squeeze(traces{k}(j,:,:))'),'r')
            
            %axis([-15 25 0 inf])
            
            subplot(3,3,k)
            
            imagesc(squeeze(traces{trialsOrdered(k)}(j,:,:))')
    
            title(titles{k})
           % xlabel('Time (s)')
   
        end
        
        saveas(h,['~/2presults/Figures/Neurons/Neuron_' num2str(neuronNumber)],'png')
        
        
        catch
            
        end
        neuronNumber = neuronNumber +1;
    
    
    end 
            
            


end
