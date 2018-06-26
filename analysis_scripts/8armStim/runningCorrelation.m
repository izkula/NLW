 clear; clc

fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat' };

f = dir('~/2pdata/20180507/');

p = {}
for i = 3:numel(f)
    filename = fullfile(f(i).folder, f(i).name)

    %load metadata
     load(fullfile(filename, 'metadata.mat'))
      counter = 1;
    for j = 1:3

    load(fullfile(filename, fnames{j}))

    %%
    C = fillmissing(C','linear')';

    %%
    for k = 1:size(C,1)
        p{i-2}(counter) = corr(C(k,:)',meta.run);
        counter = counter + 1;
    end


    end
    
    
    h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 1 1];
    hist(p{i-2})
    
    axis([-.3 .3 -inf inf])
    xlabel('Pearson Corr.')
    ylabel('Num of Neurons')
    set(gca,'XTick', [-0.3 0 .3], 'FontSize',6)
    
    saveas(h1, fullfile('~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/Running/', f(i).name),'pdf')
    
    
    
    
    
    
    
    
    
    
    
    
    
end

%%
h2 = figure; h2.Units = 'inches'; h2.Position = [1 1 2 1];
for i = 1:numel(p)
    p_mean(i) = mean(p{i});
    p_sem(i) = sem(p{i}');
end

errorb(p_mean, p_sem)
hold on; bar(p_mean)
axis([-inf inf -.2 .2])
xlabel('Mouse Number')
ylabel('Pearson Correlation')
    set(gca,'FontSize',6)

    saveas(h2, fullfile('~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/Running/','Summary'),'pdf')
