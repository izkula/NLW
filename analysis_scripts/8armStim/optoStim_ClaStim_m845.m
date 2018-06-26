
load('~/2presults/20180417/m845/metadata.mat')


    %get the trial starts
    s = diff(meta.raw.sync)';
    s = find((s>3.90)); %trial start indices
    
 nFramesBefore = 3*meta.nidaqSampleRate;
 nFramesStim = 2*meta.nidaqSampleRate;
 nFramesAfter = 3*meta.nidaqSampleRate;
 nFramesTotal = nFramesBefore + nFramesAfter + nFramesStim;
    
fp_t = zeros(10, nFramesBefore + nFramesStim + nFramesAfter);
for i = 1:numel(s)
    fp_t(i,:) = meta.fp(s(i) - nFramesBefore : s(i) + nFramesStim + nFramesAfter-1)   ;
end

%mean
fp_t = fp_t - mean(fp_t(:,1:nFramesBefore),2);

%%
h1 = figure; h1.Units = 'inches'; h1.Position = [2 2 1.5 1.5];
hold on;

fp_t_red = fp_t(1:100:end);

plot(fp_t','k','LineWidth',.2)
shadedErrorBar([],mean(fp_t),sem(fp_t))

secTick = 1; %seconds for ticks
x = (1:nFramesTotal)*1/meta.nidaqSampleRate;
x = x(1:meta.nidaqSampleRate*secTick:end);
xticks = 1:nFramesTotal;
xticks = xticks(1:meta.nidaqSampleRate*secTick:end);

set(gca,'XTick',xticks,'XTickLabel',round(x), 'FontSize', 6)

rectangle('Position',[nFramesBefore 0.08 nFramesStim 0.01],'FaceColor','r')

xlabel 'Time (s)'
ylabel 'dF (F - baseline)'
saveas(h1,'~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/PSTH_Cla_FP', 'pdf')