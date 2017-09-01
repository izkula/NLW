function run_plotLDA(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'));

dfof=dfof_sub;
%time=[0:frametime:(length(dfof)-1)*frametime];

% flag='tone';
flag='press';

% find behavioral stimulus starts
for a=1:length(tone_start)
    tone_start_2p(a)=find(time<=time_behav(tone_start(a)),1,'last');
end
for a=1:length(lever_extend)
    lever_extend_2p(a)=find(time<=time_behav(lever_extend(a)),1,'last');
end
for a=1:length(lever_retract)
    lever_retract_2p(a)=find(time<=time_behav(lever_retract(a)),1,'last');
end

% % only keep lever retractions that followed a press, not a miss
% lever_retract_2p(missed_trials)=[];

% include this data in the LDA
plot_pre=2;
plot_post=5;
pre_ind=round(plot_pre/frametime);
post_ind=round(plot_post/frametime);
total_samples=pre_ind+post_ind+1;
if plot_pre>0
    time_trial=time(1:total_samples)-time(pre_ind);
else
    time_trial=time(1:total_samples);
end

if strmatch('tone',flag,'exact')
    pressed_trials=sort([reward_trials shock_trials]);
    % delete pressed trials with reaction time >1s
    reaction_times=time_behav(lever_retract(pressed_trials))-time_behav(lever_extend(pressed_trials));
    del_long=find(reaction_times>1);
    pressed_trials(del_long)=[];
    X=[];
    for a=1:length(tone_start_2p)
        X=[X dfof(:,tone_start_2p(a)-pre_ind:tone_start_2p(a)+post_ind)];
    end
    savename=strcat(currpath,'/results/',fn,'_lda_tone.mat');
else
    X=[];
    for a=1:length(lever_retract_2p)
        X=[X dfof(:,lever_retract_2p(a)-pre_ind:lever_retract_2p(a)+post_ind)];
    end
    savename=strcat(currpath,'/results/',fn,'_lda_rewardvshock.mat');
end

% smooth and normalize X
for a=1:size(X,1)
    X(a,:)=smooth(X(a,:),15);
    maxval=max(X(a,:));
    X(a,:)=X(a,:)./maxval;
end

% find active cells
del=[];
for a=1:size(X,1)
    currcell=X(a,:);
    threshold=6*mad(currcell,1);
    
    temp=find(currcell>threshold);
    
    if isempty(temp)
        del=[del a];
    end
    
%     currcell(currcell<threshold)=0;
%     X(a,:)=currcell;
end
X(del,:)=[];
% X=X(all_cells,:);
numcells=size(X,1);

% create matrix for cov
X_all=zeros(length(lever_retract_2p),length(time_trial),size(X,1));

for c=1:size(X,1)
    currcell=X(c,:);
    currcell_trials=reshape(currcell,length(time_trial),length(lever_retract_2p));
    for a=1:length(lever_retract_2p)
        curr_trial=currcell_trials(:,a)';
        X_all(a,:,c)=curr_trial;
    end
end

% perform jackknife cross-validation
num_reward=length(reward_trials);
num_shock=length(shock_trials);
num_trials=num_reward+num_shock;

reward_traj_test=[];
shock_traj_test=[];
w_cells_all=[];
for z=1:length(lever_retract_2p)
    
    test_trial=z;
    train_reward=[];
    train_shock=[];
    
    if ismember(test_trial,shock_trials)
        test_type='shock';
        train_reward=reward_trials;
        train_shock=shock_trials(~ismember(shock_trials,test_trial));
    else
        test_type='reward';
        train_shock=shock_trials;
        train_reward=reward_trials(~ismember(reward_trials,test_trial));
    end
    
    X_reward_train=[];
    for a=1:length(train_reward)
        curr_trial=train_reward(a);
        curr_sig=squeeze(X_all(curr_trial,:,:));
        X_reward_train=[X_reward_train; curr_sig];
    end
    
    X_shock_train=[];
    for a=1:length(train_shock)
        curr_trial=train_shock(a);
        curr_sig=squeeze(X_all(curr_trial,:,:));
        X_shock_train=[X_shock_train; curr_sig];
    end

    cov_shock=cov(X_shock_train);
    cov_reward=cov(X_reward_train);
    mean_shock=mean(X_shock_train)';
    mean_reward=mean(X_reward_train)';

    w=(cov_reward+cov_shock)\(mean_reward-mean_shock);
    w_cells=[w [1:length(w)]'];
    w_cells=sortrows(w_cells);
%     w_cells(:,1)=w_cells(:,1)/abs(max(w_cells(:,1)));
    w_cells_all(:,:,z)=w_cells;
    
    if strmatch('reward',test_type)
        current_sig=squeeze(X_all(test_trial,:,:));
        reward_traj_test=[reward_traj_test; w'*current_sig'];
    else
        current_sig=squeeze(X_all(test_trial,:,:));
        shock_traj_test=[shock_traj_test; w'*current_sig'];
    end
end

figure(1);
%subplot(2,1,1);
for a=1:num_shock
%     plot(time_trial,shock_traj_test(a,:),'--r');
%     hold on
    shock_traj_test(a,:)=smooth(shock_traj_test(a,:),15);
    shock_test(a)=mean(shock_traj_test(a,:));
end
shadedErrorBar(time_trial,mean(shock_traj_test),std(shock_traj_test)/sqrt(size(shock_traj_test,1)),'r');
hold on
for a=1:num_reward
%     plot(time_trial,reward_traj_test(a,:),'--b');
%     hold on
    reward_traj_test(a,:)=smooth(reward_traj_test(a,:),15);
    reward_test(a)=mean(reward_traj_test(a,:));
end
shadedErrorBar(time_trial,mean(reward_traj_test),std(reward_traj_test)/sqrt(size(reward_traj_test,1)),'b');
plot(time_trial,zeros(1,length(time_trial)),'-k');
hold off
% subplot(2,1,2);
% figure;
% plot([1:num_shock],shock_test,'.-r');
% hold on
% plot([1:num_reward],reward_test,'.-b');
% plot([1:max([num_reward num_shock])],zeros(1,max([num_reward num_shock])),'-k');
% hold off

for a=1:num_shock+num_reward
    if a<=num_shock
        trial_types{a}='shock';
    else
        trial_types{a}='reward';
    end
end

% [Xauc,Yauc,~,AUC]=perfcurve(trial_types,[shock_test reward_test],'shock');
% 
% figure(2);
% plot(Xauc,Yauc);
% xlabel('False positive rate');
% ylabel('True positive rate');
% text(0.9,0.1,num2str(AUC),'HorizontalAlignment','c');
% hold off
% 
% % missed_cells=w_cells_all(end-4:end,2,:);
% % missed_cells=reshape(missed_cells,1,size(missed_cells,1)*size(missed_cells,3));
% % figure;
% % missed_hist=histogram(missed_cells,'binwidth',1);
% % missed_counts=[missed_hist.Values' [1:length(missed_hist.Values)]'];
% % missed_counts=sortrows(missed_counts,1);
% % top5_missed_inds=missed_counts(end-4:end,2);
% % top5_missed=missed_hist.BinEdges(top5_missed_inds);
% % 
% % pressed_cells=w_cells_all(1:5,2,:);
% % pressed_cells=reshape(pressed_cells,1,size(pressed_cells,1)*size(pressed_cells,3));
% % figure;
% % pressed_hist=histogram(pressed_cells,'binwidth',1);
% % pressed_counts=[pressed_hist.Values' [1:length(pressed_hist.Values)]'];
% % pressed_counts=sortrows(pressed_counts,1);
% % top5_pressed_inds=pressed_counts(end-4:end,2);
% % top5_pressed=pressed_hist.BinEdges(top5_pressed_inds);
% 
% p_auc=[];
% numshuff=1000;
% AUC_shuff=zeros(1,numshuff);
% for s=1:numshuff
%     shock_traj_test=[];
%     reward_traj_test=[];
% %     pressed_traj_train=[];
% %     missed_traj_train=[];
%     
%     order=randperm(num_trials);
%     shock_trials_shuff=order(1:num_shock);
%     shock_trials_shuff=sort(shock_trials_shuff);
%     reward_trials_shuff=order(num_shock+1:end);
%     reward_trials_shuff=sort(reward_trials_shuff);
%     
%     for z=1:length(lever_retract_2p)
%         
%         test_trial=z;
%         train_reward=[];
%         train_shock=[];
%         
%         if ismember(test_trial,shock_trials_shuff)
%             test_type='shock';
%             train_reward=reward_trials_shuff;
%             train_shock=shock_trials_shuff(~ismember(shock_trials_shuff,test_trial));
%         else
%             test_type='reward';
%             train_reward=reward_trials_shuff(~ismember(reward_trials_shuff,test_trial));
%             train_shock=shock_trials_shuff;
%         end
%         
%         X_reward_train=[];
%         for a=1:length(train_reward)
%             curr_trial=train_reward(a);
%             curr_sig=squeeze(X_all(curr_trial,:,:));
%             X_reward_train=[X_reward_train; curr_sig];
%         end
%         
%         X_shock_train=[];
%         for a=1:length(train_shock)
%             curr_trial=train_shock(a);
%             curr_sig=squeeze(X_all(curr_trial,:,:));
%             X_shock_train=[X_shock_train; curr_sig];
%         end
%         
%         cov_shock=cov(X_shock_train);
%         cov_reward=cov(X_reward_train);
%         mean_shock=mean(X_shock_train)';
%         mean_reward=mean(X_reward_train)';
%         
%         w=(cov_reward+cov_shock)\(mean_reward-mean_shock);
%         w_cells=[w [1:length(w)]'];
%         w_cells=sortrows(w_cells);
%          
% %         for a=1:length(train_pressed)
% %             curr_trial=pressed_trials_shuff(train_pressed(a));
% %             current_sig=squeeze(X_all(curr_trial,:,:));
% %             pressed_traj_train=[pressed_traj_train; w'*current_sig'];
% %         end
% %         
% %         for a=1:length(train_missed)
% %             curr_trial=missed_trials_shuff(train_missed(a));
% %             current_sig=squeeze(X_all(curr_trial,:,:));
% %             missed_traj_train=[missed_traj_train; w'*current_sig'];
% %         end
%         
%         if strmatch('shock',test_type)
%             current_sig=squeeze(X_all(test_trial,:,:));
%             shock_traj_test=[shock_traj_test; w'*current_sig'];
%         else
%             current_sig=squeeze(X_all(test_trial,:,:));
%             reward_traj_test=[reward_traj_test; w'*current_sig'];
%         end
%     end
%     
%     shock_test=[];
%     for a=1:num_shock
%         shock_traj_test(a,:)=smooth(shock_traj_test(a,:),15);
%         shock_test(a)=mean(shock_traj_test(a,:));
%         %search=sort(pressed_traj_test(a,:),2,'descend');
% %         pressed_test(a)=mean(pressed_traj_test(a,1:20));
%     end
%     reward_test=[];
%     for a=1:num_reward
%         reward_traj_test(a,:)=smooth(reward_traj_test(a,:),15);
%         reward_test(a)=mean(reward_traj_test(a,:));
%         %search=sort(missed_traj_test(a,:),2,'ascend');
%         %         missed_test(a)=mean(missed_traj_test(a,1:20));
%     end
%     
%     trial_types_shuff=[];
%     for a=1:num_shock+num_reward
%         if a<=num_shock
%             trial_types_shuff{a}='shock';
%         else
%             trial_types_shuff{a}='reward';
%         end
%     end
%     
%     [~,~,~,AUC_shuff(s)]=perfcurve(trial_types_shuff,[shock_test reward_test],'shock');
%     display(s);
% end
% 
% fake=length(find(AUC_shuff>=AUC));
% p_auc=fake/numshuff;
% 
% display(numcells);
% display(p_auc);
% display(AUC);
% 
% save(savename,'AUC','p_auc');