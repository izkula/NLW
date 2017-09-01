% for a=1:size(mtrials_missed_tone,1)
% figure(1);
% subplot(1,4,1);
% shadedErrorBar(time_trial,mtrials_reward_press(a,:),semtrials_reward_press(a,:));
% title(strcat('Cell ',num2str(a),'Reward'));
% subplot(1,4,2);
% shadedErrorBar(time_trial,mtrials_shock_press(a,:),semtrials_shock_press(a,:));
% title('Shock');
% subplot(1,4,3);
% shadedErrorBar(time_trial,mtrials_pressed_tone(a,:),semtrials_pressed_tone(a,:));
% title('Tone pressed');
% subplot(1,4,4);
% shadedErrorBar(time_trial,mtrials_missed_tone(a,:),semtrials_missed_tone(a,:));
% title('Tone missed');
% pause;
% end


for a=1:size(mtrials_missed_tone,1)
figure(1);
subplot(2,4,1);
shadedErrorBar(time_trial,mtrials_reward_press(a,:),semtrials_reward_press(a,:));
hold on
plot([0 0],[-0.2 0.2],'--k');
plot([1 1],[-0.2 0.2],'--k');
hold off
% ylim([-0.2 0.2])
xlim([-2 5])
title(strcat('Cell ',num2str(a),'Reward'));
subplot(2,4,2);
shadedErrorBar(time_trial,mtrials_shock_press(a,:),semtrials_shock_press(a,:));
hold on
plot([0 0],[-0.2 0.2],'--k');
plot([1 1],[-0.2 0.2],'--k');
plot([2 2],[-0.2 0.2],'--k');
hold off
% ylim([-0.2 0.2])
xlim([-2 5])
title('Shock');
subplot(2,4,3);
shadedErrorBar(time_trial,mtrials_pressed_tone(a,:),semtrials_pressed_tone(a,:));
hold on
plot([0 0],[-0.2 0.2],'--k');
plot([5 5],[-0.2 0.2],'--k');
hold off
% ylim([-0.2 0.2])
xlim([-2 5])
title('Tone pressed');
subplot(2,4,4);
shadedErrorBar(time_trial,mtrials_missed_tone(a,:),semtrials_missed_tone(a,:));
hold on
plot([0 0],[-0.2 0.2],'--k');
plot([5 5],[-0.2 0.2],'--k');
hold off
% ylim([-0.2 0.2])
xlim([-2 5])
title('Tone missed');

subplot(2,4,5);
for b=1:size(dfof_trials_reward_press,2)
    plot(time_trial,squeeze(dfof_trials_reward_press(a,b,:)));
%     ylim([-0.5 1])
    xlim([-2 5])
    hold on
end
plot([0 0],[-0.5 1],'--k');
plot([1 1],[-0.5 1],'--k');
hold off
subplot(2,4,6);
for b=1:size(dfof_trials_shock_press,2)
    plot(time_trial,squeeze(dfof_trials_shock_press(a,b,:)));
%     ylim([-0.5 1])
    xlim([-2 5])
    hold on
end
plot([0 0],[-0.5 1],'--k');
plot([1 1],[-0.5 1],'--k');
plot([2 2],[-0.5 1],'--k');
hold off
subplot(2,4,7);
for b=1:size(dfof_trials_pressed_tone,2)
    plot(time_trial,squeeze(dfof_trials_pressed_tone(a,b,:)));
%     ylim([-0.5 1])
    xlim([-2 5])
    hold on
end
plot([0 0],[-0.5 1],'--k');
plot([5 5],[-0.5 1],'--k');
hold off
subplot(2,4,8);
for b=1:size(dfof_trials_missed_tone,2)
    plot(time_trial,squeeze(dfof_trials_missed_tone(a,b,:)));
%     ylim([-0.5 1])
    xlim([-2 5])
    hold on
end
plot([0 0],[-0.5 1],'--k');
plot([5 5],[-0.5 1],'--k');
hold off
pause;
end