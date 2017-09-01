close all
clear all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
filedate='20160429';
fn='GRIN_pfcnac5_tonerewardshock_day5';

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/results/',fn,'.mat'));
load(strcat(currpath,'/results/',fn,'_lda_press_allcells.mat'));



figure;
hold on;
for a=1:5
    subplot(2,5,a);
    for b=1:size(dfof_trials_pressed_tone,2)
        plot(time_trial,squeeze(dfof_trials_pressed_tone(top5_pressed(a),b,:)),'b');
        hold on
    end
    for b=1:size(dfof_trials_missed_tone,2)
        plot(time_trial,squeeze(dfof_trials_missed_tone(top5_pressed(a),b,:)),'r');
        hold on
    end
%     shadedErrorBar(time_trial,mtrials_pressed_tone(top5_pressed(a),:),semtrials_pressed_tone(top5_pressed(a),:),'b');
%     hold on
%     shadedErrorBar(time_trial,mtrials_missed_tone(top5_pressed(a),:),semtrials_missed_tone(top5_pressed(a),:),'r');
    xlim([-2 10])


    subplot(2,5,a+5);
    for b=1:size(dfof_trials_reward_press,2)
        plot(time_trial,squeeze(dfof_trials_reward_press(top5_pressed(a),b,:)),'b');
        hold on
    end
    for b=1:size(dfof_trials_shock_press,2)
        plot(time_trial,squeeze(dfof_trials_shock_press(top5_pressed(a),b,:)),'r');
        hold on
    end
   
%     shadedErrorBar(time_trial,mtrials_reward_press(top5_pressed(a),:),semtrials_reward_press(top5_pressed(a),:),'b');
%     hold on
%     shadedErrorBar(time_trial,mtrials_shock_press(top5_pressed(a),:),semtrials_shock_press(top5_pressed(a),:),'r');
    xlim([-2 10])
end
title('Pressed Cells');


figure;
hold on;
for a=1:5
    subplot(2,5,a);
%     for b=1:size(dfof_trials_pressed_tone,2)
%         plot(time_trial,squeeze(dfof_trials_pressed_tone(top5_pressed(a),b,:)),'b');
%         hold on
%     end
%     for b=1:size(dfof_trials_missed_tone,2)
%         plot(time_trial,squeeze(dfof_trials_missed_tone(top5_pressed(a),b,:)),'r');
%         hold on
%     end
    shadedErrorBar(time_trial,mtrials_pressed_tone(top5_pressed(a),:),semtrials_pressed_tone(top5_pressed(a),:),'b');
    hold on
    shadedErrorBar(time_trial,mtrials_missed_tone(top5_pressed(a),:),semtrials_missed_tone(top5_pressed(a),:),'r');
    xlim([-2 10])


    subplot(2,5,a+5);
%     for b=1:size(dfof_trials_reward_press,2)
%         plot(time_trial,squeeze(dfof_trials_reward_press(top5_pressed(a),b,:)),'b');
%         hold on
%     end
%     for b=1:size(dfof_trials_shock_press,2)
%         plot(time_trial,squeeze(dfof_trials_shock_press(top5_pressed(a),b,:)),'r');
%         hold on
%     end
   
    shadedErrorBar(time_trial,mtrials_reward_press(top5_pressed(a),:),semtrials_reward_press(top5_pressed(a),:),'b');
    hold on
    shadedErrorBar(time_trial,mtrials_shock_press(top5_pressed(a),:),semtrials_shock_press(top5_pressed(a),:),'r');
    xlim([-2 10])
end
title('Pressed Cells');





figure;
hold on;
for a=1:5
    subplot(2,5,a);
%     for b=1:size(dfof_trials_pressed_tone,2)
%         plot(time_trial,squeeze(dfof_trials_pressed_tone(top5_missed(a),b,:)),'b');
%         hold on
%     end
%     for b=1:size(dfof_trials_missed_tone,2)
%         plot(time_trial,squeeze(dfof_trials_missed_tone(top5_missed(a),b,:)),'r');
%         hold on
%     end
    
    shadedErrorBar(time_trial,mtrials_pressed_tone(top5_missed(a),:),semtrials_pressed_tone(top5_missed(a),:),'b');
    hold on
    shadedErrorBar(time_trial,mtrials_missed_tone(top5_missed(a),:),semtrials_missed_tone(top5_missed(a),:),'r');

    xlim([-2 10])

    subplot(2,5,a+5);
%     for b=1:size(dfof_trials_reward_press,2)
%         plot(time_trial,squeeze(dfof_trials_reward_press(top5_missed(a),b,:)),'b');
%         hold on
%     end
%     for b=1:size(dfof_trials_shock_press,2)
%         plot(time_trial,squeeze(dfof_trials_shock_press(top5_missed(a),b,:)),'r');
%         hold on
%     end

    shadedErrorBar(time_trial,mtrials_reward_press(top5_missed(a),:),semtrials_reward_press(top5_missed(a),:),'b');
    hold on
    shadedErrorBar(time_trial,mtrials_shock_press(top5_missed(a),:),semtrials_shock_press(top5_missed(a),:),'r');

    xlim([-2 10])
end
title('Missed Cells');





figure;
hold on;
for a=1:5
    subplot(2,5,a);
    for b=1:size(dfof_trials_pressed_tone,2)
        plot(time_trial,squeeze(dfof_trials_pressed_tone(top5_missed(a),b,:)),'b');
        hold on
    end
    for b=1:size(dfof_trials_missed_tone,2)
        plot(time_trial,squeeze(dfof_trials_missed_tone(top5_missed(a),b,:)),'r');
        hold on
    end
    
%     shadedErrorBar(time_trial,mtrials_pressed_tone(top5_missed(a),:),semtrials_pressed_tone(top5_missed(a),:),'b');
%     hold on
%     shadedErrorBar(time_trial,mtrials_missed_tone(top5_missed(a),:),semtrials_missed_tone(top5_missed(a),:),'r');

    xlim([-2 10])

    subplot(2,5,a+5);
    for b=1:size(dfof_trials_reward_press,2)
        plot(time_trial,squeeze(dfof_trials_reward_press(top5_missed(a),b,:)),'b');
        hold on
    end
    for b=1:size(dfof_trials_shock_press,2)
        plot(time_trial,squeeze(dfof_trials_shock_press(top5_missed(a),b,:)),'r');
        hold on
    end

%     shadedErrorBar(time_trial,mtrials_reward_press(top5_missed(a),:),semtrials_reward_press(top5_missed(a),:),'b');
%     hold on
%     shadedErrorBar(time_trial,mtrials_shock_press(top5_missed(a),:),semtrials_shock_press(top5_missed(a),:),'r');

    xlim([-2 10])
end
title('Missed Cells');