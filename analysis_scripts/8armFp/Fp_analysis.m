clear; clc
init_samz

fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8armFp/Figures/'
f = dir('~/2pdata/20180522/')


%An IMPORTANT note on the trial structure%
%Basically, the frames are unreliable to indicate when a trial starts. This
%is because the roration of the 8arm pushes the trial length over the 30 s
%allowable trial length. Thus,even though the stimulation period begins
%correctly, the object is stil moving towards the mouse during the
%'stationary' period. A reliable way to align trials is to use the
%meta.framesAlined(2*trialnumber) to get teh frame when both the (1) stim e
%nds and (2) retreat begins. I have measured many trial values for several
%mice and trials, and after trial 1, the stim is consistently ~310 frames.
%The stationary lenght is consistently ~212 frames. The trial start to
%trial start is consistently 926ish frames. 

%So, I am going to use the end of the stim period as the algining frame
%number. 

nFramesStim = round(10*30.98)
nFramesBeforeStim = 3*31
nFramesBeforeStimEnd =  nFramesStim + nFramesBeforeStim;
nFramesAfterStimEnd = 3*31;
nFramesTrial = nFramesBeforeStimEnd + nFramesAfterStimEnd


FPT = []; RUNT = []; M = {} 
for i = 3:numel(f)
    foldername = fullfile('~/2pdata/20180522/', f(i).name);
    if contains(foldername,'fp') 
    load(fullfile(foldername, 'metadata.mat')) %load metadata
    counter = 1;
    
    fp_t = nan(1,nFramesTrial , size(meta.frames,1)/2);
    run_t = nan(1,nFramesTrial, size(meta.frames,1)/2);
    
    %zscore
    fp_z = zscore(meta.fp);
    
    for j = 1:(size(meta.framesAligned,1)/2)  %for each trial      
        try
         stimEndFrame =  meta.framesAligned(2*j) %define trial marker
     
         fp_t(:,:,j) = fp_z(stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1,1);
         run_t(:,:,j) = meta.run(stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1,1);
         
         % eye_t(:,:,j) = meta.eye(stimEndFrame - nFramesBeforeStimEnd : stimEndFrame + nFramesAfterStimEnd-1);
         
        catch
            fprintf('messed up' )
            i
            j
                
        end      
    end
    FPT = [FPT;fp_t]
    RUNT = [RUNT; run_t];
    M{i} = meta; 
    end
end

    

clearvars -except FPT RUNT M nFramesBeforeStim nFramesStim nFramesBeforeStimEnd nFramesAfterStimEnd nFramesTrial

save('~/2presults/8armFP/allMiceTrial.mat', '-v7.3')

%% Housekeeping and Loading in things

init_samz
%fnames = {'AVG_z0_neuron.mat' 'AVG_z1_neuron.mat' 'AVG_z2_neuron.mat'}
fig_dir = '/home/svesuna/Dropbox/2p_Claustrum_Shared/2p/Results/8armFP/Figures/'
%load('~/2presults/8armFP/allMiceTrial.mat')

%hardcode in the trials
stimTrials = repmat([ones(1,8) zeros(1,8)],1,4)
noStimTrials =abs(stimTrials -1);
objectTrials = repmat([1 0 1 0 0 1 0 0],1,8)
socialTrials = repmat([0 1 0 0 1 0 1 0],1,8)
blankTrials = repmat([0 0 0 1 0 0 0 0],1,8)
toyTrials = repmat([0 0 0 0 0 0 0 1], 1,8)

allTrials = [objectTrials.*noStimTrials; objectTrials.*stimTrials;...
    socialTrials.*noStimTrials; socialTrials.*stimTrials;...
    toyTrials.*noStimTrials; toyTrials.*stimTrials;...
    blankTrials.*noStimTrials; blankTrials.*stimTrials];

%% Calculate the average FP number fo each trial type for each mouse

%FPT2 = FPT - mean(FPT(:,1:nFramesBeforeStim,:),2)

FPT_avg = [];
for i = 1:4*2 %types of trials we want to include
    x1 = FPT(:,:,find(allTrials(i,:)));
 
    
    x2 = nanmean(x1(:,nFramesBeforeStim:nFramesBeforeStim +nFramesStim,:),2);
  
    x3 = nanmean(x2,3);
    
    FPT_avg(:,i) = x3;
    
end

%% statistics
[P ANOVATAB STATS] = anova1(FPT_avg);
C = multcompare(STATS);

%% plot
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 3 3];
errorbar(mean(FPT_avg,1),sem(FPT_avg));
hold on;
bar(mean(FPT_avg));

for i = 1:8
    scatter(i*ones(7,1),FPT_avg(:,i), 2, 'MarkerFaceColor','k','MarkerEdgeColor','k')
end
    
set(gca,'XTick',[1:8],'XTickLabel',{'Object' 'Object Stim' 'Social' 'Social Stim' 'Toy Mouse' 'Toy Mouse Stim' 'Blank' 'Blank Stim'},...
    'FontSize',6,'XTickLabelRotation', 45)
ylabel( 'Photometry (zscore)')

saveas(h1,fullfile(fig_dir, 'barplot'),'pdf')
    
    















