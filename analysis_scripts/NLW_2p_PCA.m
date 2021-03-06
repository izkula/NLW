clear; close all; clc
init_samz
%%The script imports trialtraces and then plots all neurons%%

%Logic
doIndividualStimulusFigures = 1

f = dir('~/2presults/8arm/')

%aggregates trace data from all mice and trials and creates one trial
%averaged matrix


for i = 3:numel(f)
    dirName = fullfile(f(i).folder,f(i).name);
    
    %load neuron trial trace data
    load(fullfile(dirName, 'neuron.mat'));
    load(fullfile(dirName, 'metadata.mat'));

    
    
    [COEFF, SCORE, LATENT, TSQ, EXPLAINED] = pca(C');
    
    
    
    plot(SCORE(:,1:3))
    hold on
    PlotVerticalLines(meta.frames)
end
