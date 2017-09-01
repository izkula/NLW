function press_vs_miss_barplot(pathname,filedate,fn)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/results/',fn,'.mat'));
load(strcat(currpath,'/behavior/',fn,'.mat'));
load(strcat(currpath,'/regression/',fn,'_shockonly.mat'));
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
dfof=dfof_sub;

behav_cells=shock_only;

t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial<5,1,'last');

trial_type3=dfof_trials_missed_tone;
numtrials3=size(trial_type3,2);

% only look at trials where he pressed within 1s
pressed_trials=sort([reward_trials shock_trials]);
reaction_times=time_behav(lever_retract(pressed_trials))-time_behav(lever_extend(pressed_trials));
figure; plot([1:length(reaction_times)],reaction_times,'o');
del_long=find(reaction_times>1);
hold on
plot(del_long,reaction_times(del_long),'x');
shock_inds=ismember(pressed_trials,shock_trials);
shock_inds=find(shock_inds==1);
plot(shock_inds,reaction_times(shock_inds),'og');

trial_type4=dfof_trials_pressed_tone;
trial_type4(:,del_long,:)=[];
numtrials4=size(trial_type4,2);

stdthresh=2;

for a=1:length(behav_cells)
   currdfof=squeeze(trial_type3(behav_cells(a),:,:));
   active_trials=[];
   for b=1:numtrials3
       currtrial=currdfof(b,t_0:t_5);
       search=find(currtrial>stdthresh);
       if ~isempty(search)
           active_trials=[active_trials b];
       end
   end
   peaks_missed(a)=length(active_trials)/numtrials3;
end

for a=1:length(behav_cells)
   currdfof=squeeze(trial_type4(behav_cells(a),:,:));
   if size(trial_type4,2)<2
       currdfof=currdfof';
   end
   active_trials=[];
   for b=1:numtrials4
       currtrial=currdfof(b,t_0:t_5);
       search=find(currtrial>stdthresh);
       if ~isempty(search)
           active_trials=[active_trials b];
       end
   end
   peaks_pressed(a)=length(active_trials)/numtrials4;
end

p_val=signrank(peaks_missed,peaks_missed);

display(p_val);
m_missed=mean(peaks_missed);
m_pressed=mean(peaks_pressed);
si=(peaks_missed-peaks_pressed)./(peaks_missed+peaks_pressed);
m_si=mean(si(~isnan(si)));
display(m_missed);
display(m_pressed);
display(m_si);

figure;
bar(1,mean(peaks_missed),'r');
hold on
bar(2,mean(peaks_pressed),'b');
legend('Missed','Pressed');
errorbar(1,mean(peaks_missed),std(peaks_missed)/sqrt(length(peaks_missed)),'r');
errorbar(2,mean(peaks_pressed),std(peaks_pressed)/sqrt(length(peaks_pressed)),'b');

savename=strcat(currpath,'/regression/',fn,'_peaks.mat');
save(savename,'peaks_missed','peaks_pressed','p_val');

end
