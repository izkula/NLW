clear; close all; clc

init_samz

%load data
try
    z0 = load('~/2presults/Running/m876_running_003/z0_neuron.mat');
    z1 = load('~/2presults/Running/m876_running_003/z1_neuron.mat');
    z2 = load('~/2presults/Running/m876_running_003/z2_neuron.mat');
    load('~/2presults/Running/m876_running_003/metadata.mat')
    fprintf('All three planes and meta data loaded')
catch
    error('Couldnt load data or meta')
end

%%
%Aggregate Neurons for Plotting
C = [z0.C ; z1.C; z2.C];
S = [z0.S ; z1.S; z2.S];

%Crop Both
C = C(:,meta.frameStart:meta.frameEnd-1);
S = S(:,meta.frameStart:meta.frameEnd-1);

%Scale Each Neuron
C_z = zscore(C')';

%Interpolate Movement Length to Make N Frames
R = meta.run;

%%
if isempty(find(isnan(C)))
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(C);
else
   C1 = fillmissing(C,'previous',2);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(C1);
end

figure; hold on
scatter(R,COEFF(:,1))
scatter(R,COEFF(:,2))
xlabel 'Movement'
ylabel 'PCA coeff'


%% CORRELATION BETWEEN NEURON AND MOVEMENT%
%Calculate correlation between each neuron and movement
corr_CR = [];

for i = 1:size(C_z,1)
    corr_CR(i) = corr(C_z(i,:)', R) ;
end

[corr_CR_sorted, I] = sort(corr_CR);

%Reorder neurons 
C_sorted = sortrows(C, I,'descend');
S_sorted = sortrows(S,I,'descend') ;
corr_CR_sorted = sortrows(corr_CR(I));



%% PLOT
width = 4; height = 4;

%Neuron Figure
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 width height];
subplot(6,6,[7:11, 13:17, 19:23, 25:29, 31:35])
imagesc(C_sorted)

ylabel(' Claustrum Neuron Number')
%Make X Axis Correctly
xlabel('Time (s)')
x = (1:length(C))*(1/30.98);
x = x(1:1860:end)
xticks = 1:length(C)
xticks = xticks(1:1860:end)
set(gca,'XTick', xticks, 'XTickLabel',round(x),'XTickLabelRotation', 45,...
    'FontSize',5)

%Movement Figure
subplot(6,6,1:5)
plot(meta.runCountsSmooth)
axis([-inf inf 0 25])
set(gca,'XTick',[],'FontSize',5)
ylabel('Movement au')

%Correlation Each Neuron
subplot(6,6,[12 18 24 30 36]);
barh(corr_CR_sorted)
xlabel 'Correlation'
set(gca,'XTick',[-1 0 1],'YTick',[],'FontSize',5)

%Correlation Histogram
subplot(6,6,6)
hist(corr_CR)
set(gca,'FontSize',5)
ylabel 'n Neurons'

saveas(h1, '~/Dropbox/2p_Claustrum_Shared/2p/Results/Running/Figures/MovementRaw','pdf')


%%




clear; close all; clc

init_samz

%load data
try
    z0 = load('~/2presults/Running/m876_running_003/z0_neuron.mat');
    z1 = load('~/2presults/Running/m876_running_003/z1_neuron.mat');
    z2 = load('~/2presults/Running/m876_running_003/z2_neuron.mat');
    load('~/2presults/Running/m876_running_003/metadata.mat')
    fprintf('All three planes and meta data loaded')
catch
    error('Couldnt load data or meta')
end


%%
%
hist(corr_CR)
xlabel 'Pearson Correlation with Movement'
ylabel 'Number of Neurons'