function regression_barplot(pathname,filedate,fn)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/results/',fn,'.mat'));
load(strcat(currpath,'/behavior/',fn,'.mat'));
load(strcat(currpath,'/regression/',fn,'_shockonly.mat'));

% trial_type=dfof_trials_shock_press;
% 
% numtrials=size(trial_type,2);
% numtimepoints=length(time_trial);

% t_1=find(time_trial>=1,1,'first');
% t_3=find(time_trial>=3,1,'first');
% behav=zeros(1,length(time_trial));
% behav(t_1:t_3)=5;
% % behav_regress=repmat(behav,1,numtrials);
% 
% % convolve behavior regressor with exponential
% tau=0.25;
% cirf=exp(-time_trial/tau);
% cirf=cirf/sum(cirf);
% behav_conv=conv(behav,cirf);
% behav_conv=behav_conv(1:length(cirf));
% behav_regress=repmat(behav_conv,1,numtrials);
% 
% cells_type=[];
% for a=1:numcells
%     cells_type(a,:)=reshape(squeeze(trial_type(a,:,:))',1,numtimepoints*numtrials);
%     mcells_type(a,:)=mean(squeeze(trial_type(a,:,:)));
%     semcells_type(a,:)=std(squeeze(trial_type(a,:,:)))/sqrt(numtrials);
% end

behav_cells=shock_only;

% for a=1:length(behav_cells)
%     figure(111);
%     plot(cells_type(behav_cells(a),:));
%     hold on;
% end
% plot([1:numtrials*numtimepoints],behav_regress,'r','linewidth',2);
%plot(time_trial,behav);
% display(length(behav_cells));
% display(length(behav_cells)/numcells);

% trial_type2=dfof_trials_reward_press;
% numtrials2=size(trial_type2,2);
% cells_type2=[];
% for a=1:numcells
%     cells_type2(a,:)=reshape(squeeze(trial_type2(a,:,:))',1,numtimepoints*numtrials2);
%     mcells_type2(a,:)=mean(squeeze(trial_type2(a,:,:)));
%     semcells_type2(a,:)=std(squeeze(trial_type2(a,:,:)))/sqrt(numtrials2);
% end
% behav_regress2=repmat(behav_conv,1,numtrials2);
% 
% [rho2,p2]=corr(behav_regress2',cells_type2');
% overlap=find(p2(behav_cells)<0.05);


t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial>=5,1,'first');
% behav2=zeros(1,length(time_trial));
% behav2(t_0:t_5)=5;

numtimepoints=length(time_trial);

trial_type3=dfof_trials_missed_tone;
for a=1:size(trial_type3,1)
    currcell=squeeze(trial_type3(a,:,:));
    for b=1:size(currcell,1)
        currcell(b,:)=smooth(currcell(b,:),10);
    end
    trial_type3(a,:,:)=currcell;
end


numcells=size(trial_type3,1);
numtrials3=size(trial_type3,2);
cells_type3=[];
for a=1:numcells
    cells_type3(a,:)=reshape(squeeze(trial_type3(a,:,:))',1,numtimepoints*numtrials3);
    mcells_type3(a,:)=mean(squeeze(trial_type3(a,:,:)));
    semcells_type3(a,:)=std(squeeze(trial_type3(a,:,:)))/sqrt(numtrials3);
end


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

for a=1:size(trial_type4,1)
    currcell=squeeze(trial_type4(a,:,:));
    for b=1:size(currcell,1)
        currcell(b,:)=smooth(currcell(b,:),10);
    end
    trial_type4(a,:,:)=currcell;
end

trial_type4(:,del_long,:)=[];
numtrials4=size(trial_type4,2);
cells_type4=[];
for a=1:numcells
    cells_type4(a,:)=reshape(squeeze(trial_type4(a,:,:))',1,numtimepoints*numtrials4);
    mcells_type4(a,:)=mean(squeeze(trial_type4(a,:,:)));
    semcells_type4(a,:)=std(squeeze(trial_type4(a,:,:)))/sqrt(numtrials4);
end

%peaks_pressed=max(trial_type4(behav_cells,:,t_0:t_5),[],3);
peaks_pressed=mean(trial_type4(behav_cells,:,t_0:t_5),3);
peaks_pressed=mean(peaks_pressed,2);

%peaks_missed=max(trial_type3(behav_cells,:,t_0:t_5),[],3);
peaks_missed=mean(trial_type3(behav_cells,:,t_0:t_5),3);
peaks_missed=mean(peaks_missed,2);

p_val=signrank(peaks_missed,peaks_missed);

display(p_val);
m_missed=mean(peaks_missed);
m_pressed=mean(peaks_pressed);
display(m_missed);
display(m_pressed);

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
