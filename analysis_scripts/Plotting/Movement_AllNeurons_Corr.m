clear; close all; clc


%load data
try
    z0 = load('~/Dropbox/2p_Claustrum_Shared/2p/Results/Running/m876_run/z0_neuron.mat');
    z1 = load('~/Dropbox/2p_Claustrum_Shared/2p/Results/Running/m876_run/z1_neuron.mat');
    z2 = load('~/Dropbox/2p_Claustrum_Shared/2p/Results/Running/m876_run/z2_neuron.mat');
    load('~/Dropbox/2p_Claustrum_Shared/2p/Results/Running/m876_run/metadata.mat')
    fprintf('All three planes and meta data loaded')
catch
    error('Couldnt load data or meta')
end

%%
%Aggregate Neurons for Plotting
C = [z0.C ; z1.C; z2.C];
S = [z0.S ; z1.S; z2.S];

%Scale Each Neuron
C_z = zscore(C,1);

%Interpolate Movement Length to Make N Frames
R = interp1(meta.runCountsSmooth,1:length(C));

%% CORRELATION BETWEEN NEURON AND MOVEMENT%
%Calculate correlation between each neuron and movement
corr_CR = [];

for i = 1:size(C_z,1)
    corr_CR(i) = corr(C_z(i,:)', R','rows','complete') ;
end

[corr_CR_sorted, I] = sort(corr_CR);

%Reorder neurons 
C_sorted = sortrows(C, I);

%

%% PLOT
width = 3.5; height = 2;

%Neuron Figure
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 width height];
subplot(6,1,[2:6])
imagesc(C_sorted)

ylabel('Neuron Number')
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:length(C))*(1/30.98);
x = x(1:1860:end)
xticks = 1:length(C)
xticks = xticks(1:1860:end)
set(gca,'XTick', xticks, 'XTickLabel',round(x),'XTickLabelRotation', 45)

%Movement Figure
subplot(6,1,1)
plot(meta.runCountsSmooth)
axis([-inf inf 0 25])
set(gca,'XTick',[])
ylabel('Movement au')

%%
