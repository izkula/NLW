function press_vs_miss_barplot(pathname,filedate,fn)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/results/',fn,'.mat'));
load(strcat(currpath,'/behavior/',fn,'.mat'));
load(strcat(currpath,'/regression/',fn,'_rewardonly.mat'));
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
dfof=dfof_sub;

behav_cells=reward_only;

t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial<=5,1,'last');

trial_type3=dfof_trials_missed_tone;
numtrials3=size(trial_type3,2);

% only look at trials where he pressed within 1s
pressed_trials=sort([reward_trials shock_trials]);
reaction_times=time_behav(lever_retract(pressed_trials))-time_behav(lever_extend(pressed_trials));
figure(1); plot([1:length(reaction_times)],reaction_times,'o');
del_long=find(reaction_times>1);
hold on
plot(del_long,reaction_times(del_long),'x');
shock_inds=ismember(pressed_trials,shock_trials);
shock_inds=find(shock_inds==1);
plot(shock_inds,reaction_times(shock_inds),'og');
hold off

trial_type4=dfof_trials_pressed_tone;
trial_type4(:,del_long,:)=[];
numtrials4=size(trial_type4,2);

for a=1:length(behav_cells)
   currdfof=squeeze(trial_type3(behav_cells(a),:,:));
   if size(trial_type3(behav_cells(a),:,:),2)==1
        currdfof=currdfof';
   end
   peaks_missed(a)=mean(mean(currdfof(:,t_0:t_5),2));
end

for a=1:length(behav_cells)
   currdfof=squeeze(trial_type4(behav_cells(a),:,:));
   if size(trial_type4(behav_cells(a),:,:),2)==1
        currdfof=currdfof';
   end
   peaks_pressed(a)=mean(mean(currdfof(:,t_0:t_5),2));
end

% figure;
% for a=1:length(behav_cells)
%     plot([1 2],[peaks_missed(a) peaks_pressed(a)],'.-')
%     hold on
% end

p_val=signrank(peaks_missed,peaks_missed);


savename=strcat(currpath,'/regression/',fn,'_peaks_reward.mat');
save(savename,'peaks_missed','peaks_pressed','p_val');

end
