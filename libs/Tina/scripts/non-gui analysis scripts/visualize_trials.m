close all
clear all

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};


filepath='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
filedate='20160528';
filename='GRIN_pfcnac3_tonerewardshock_day4';

load(strcat(filepath,filedate,'/results/',filename,'.mat'));
load(strcat(filepath,filedate,'/regression/',filename,'_shockonly.mat'));

behav_cells=shock_only;

for a=1:length(behav_cells)
    for b=1:size(dfof_trials_missed_tone,2)
        dfof_trials_missed_tone(behav_cells(a),b,:)=smooth(squeeze(dfof_trials_missed_tone(behav_cells(a),b,:)),4);
    end
    for b=1:size(dfof_trials_pressed_tone,2)
        dfof_trials_pressed_tone(behav_cells(a),b,:)=smooth(squeeze(dfof_trials_pressed_tone(behav_cells(a),b,:)),4);        
    end
    for b=1:size(dfof_trials_shock_press,2)
        dfof_trials_shock_press(behav_cells(a),b,:)=smooth(squeeze(dfof_trials_shock_press(behav_cells(a),b,:)),4);
    end
    for b=1:size(dfof_trials_reward_press,2)
        dfof_trials_reward_press(behav_cells(a),b,:)=smooth(squeeze(dfof_trials_reward_press(behav_cells(a),b,:)),4);        
    end    
    
end

for a=1:length(behav_cells)
    figure(1);
%     subplot(2,1,1);
    plot(time_trial,squeeze(dfof_trials_missed_tone(behav_cells(a),:,:)),'r');
%     subplot(4,1,3);
%     plot(time_trial,mean(squeeze(dfof_trials_missed_tone(behav_cells(a),:,:))),'r');
%     hold on
%     plot(time_trial,mean(squeeze(dfof_trials_pressed_tone(behav_cells(a),:,:))),'b');
%     hold off
%     subplot(2,1,2);
    box off
    xlim([0 5])
    hold on
    plot(time_trial,squeeze(dfof_trials_pressed_tone(behav_cells(a),:,:)),'b');
%     subplot(4,1,4);
%     plot(time_trial,mean(squeeze(dfof_trials_missed_tone(behav_cells(a),:,:)))-mean(squeeze(dfof_trials_pressed_tone(behav_cells(a),:,:))),'k');
    hold off
    title(strcat('Cell:',num2str(behav_cells(a))));
    
    figure(2);
    shadedErrorBar(time_trial,smooth(mtrials_shock_press(behav_cells(a),:),4),smooth(semtrials_shock_press(behav_cells(a),:),4),'r');
    hold on
    shadedErrorBar(time_trial,smooth(mtrials_reward_press(behav_cells(a),:),4),smooth(semtrials_reward_press(behav_cells(a),:),4),'b');
    hold off
    box off
    xlim([-1 5])
    pause;
    
end