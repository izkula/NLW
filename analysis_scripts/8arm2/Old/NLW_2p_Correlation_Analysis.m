clear; close all; clc
%init_samz

init_samMac
%The script imports trialtraces and then plots all neurons%%


f = dir(fullfile(resultsPath, '8arm'))

%aggregates trace data from all mice and trials and creates one trial
%averaged matrix

trialsOrdered = [1 2 8 3 5 4 6 7]
titles = {'Obj 1' 'Soc 1' 'Empty' 'Obj 2'  'Soc 2' 'ToySoc' 'Obj3' 'Soc 3'} 

for i = 4:numel(f)
    dirName = fullfile(f(i).folder,f(i).name);
    
    %load neuron trial trace data
    load(fullfile(dirName, 'trialtraces.mat'))
    
    
    %calculate distance between neurons
    
    
    
    
    
    
    
    
    
    
    %calculate a correlation for all the data, then for each part of the
    %epoch
    
    traces_extended = reshape(traces_preliminary, nNeurons, [],1);
    
    
    
    
    
    
%     
%     tic
%     [acor, lag] = xcorr(traces_extended',0,'coeff');
%     toc
%     
%     figure
%     hist(acor,10)
%     title(f(i).name)
    
    
    
    
    
    
    
    
    
    
    
    
    
end
