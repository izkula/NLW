function run_plotPCA_v2(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'));

dfof=dfof_sub;
%time=[0:frametime:(length(dfof)-1)*frametime];

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

% create dfof for PCA
X=dfof;
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
[F,pcs,eigs,var,pvar]=pca(X);
F=F';
figure(500);
plot(pcs(:,1));
hold on
plot(pcs(:,2),'r');
plot(pcs(:,3),'g');
legend('PC1','PC2','PC3');
figure(501);
plot(pvar);

% plot this data in trajectories
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

% plot trajectories aligned to the press
% [mpca_reward_tone, sempca_reward_tone, pca_trials_reward_tone]=get_traj(F,tone_start_2p,reward_trials,'reward',time_trial,pre_ind,post_ind);
% [mpca_shock_tone, sempca_shock_tone, pca_trials_shock_tone]=get_traj(F,tone_start_2p,shock_trials,'shock',time_trial,pre_ind,post_ind);
% [mpca_miss_tone, sempca_miss_tone, pca_trials_miss_tone]=get_traj(F,tone_start_2p,missed_trials,'miss',time_trial,pre_ind,post_ind);

% get the reward trials closest to the shock trials
search_reward=reward_trials;
for a=1:length(shock_trials)
    if shock_trials(a)<search_reward(1)
        reward_nearshock(a)=search_reward(find(search_reward>shock_trials(a),1,'first'));
        search_reward(find(search_reward>shock_trials(a),1,'first'))=[];
    else
        reward_nearshock(a)=search_reward(find(search_reward<shock_trials(a),1,'last'));
        search_reward(find(search_reward<shock_trials(a),1,'last'))=[];
    end
end

shock_trials
reward_nearshock

% for a=1:length(shock_trials)
%     reward_preshock(a)=reward_trials(find(reward_trials<shock_trials(a),1,'last'));
% end

[mpca_reward, sempca_reward, pca_trials_reward]=get_traj(F,lever_retract_2p,reward_nearshock,'reward',time_trial,pre_ind,post_ind);
[mpca_shock, sempca_shock, pca_trials_shock]=get_traj(F,lever_retract_2p,shock_trials,'shock',time_trial,pre_ind,post_ind);


if ~exist(strcat(currpath,'/results/'),'dir')
    mkdir(strcat(currpath,'/results/'));
end

savename=strcat(currpath,'/results/',fn,'_pca.mat');
%save(savename,'time_trial','mpca_shock_tone','sempca_shock_tone','pca_trials_shock_tone','mpca_reward_tone','sempca_reward_tone','pca_trials_reward_tone','mpca_miss_tone','sempca_miss_tone','pca_trials_miss_tone');
save(savename,'time_trial','mpca_shock','sempca_shock','pca_trials_shock','mpca_reward','sempca_reward','pca_trials_reward');

end


function [mtrials,semtrials,F_trials]=get_traj(F,alignment,trial_type,label,time_trial,pre_ind,post_ind)

time_zero=find(time_trial>=1,1,'first');
time_five=find(time_trial>=5,1,'first');

F_trials=zeros(3,length(trial_type),length(time_trial));
for a=1:3
    for b=1:length(trial_type)
        curr_trial=trial_type(b);
        F_trials(a,b,:)=F(a,alignment(curr_trial)-pre_ind:alignment(curr_trial)+post_ind);
        F_trials(a,b,:)=smooth(F_trials(a,b,:),15);
    end
    mtrials(a,:)=mean(squeeze(F_trials(a,:,:)));
    semtrials(a,:)=std(squeeze(F_trials(a,:,:)))/sqrt(size(F_trials,2));
end

figure(123);
hold on
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

if strmatch('reward',label)
    figure(123);
    for a=1:length(trial_type)
        plot3(squeeze(F_trials(1,a,:)),squeeze(F_trials(2,a,:)),squeeze(F_trials(3,a,:)),'-g');
        hold on
                plot3(squeeze(F_trials(1,a,1)),squeeze(F_trials(2,a,1)),squeeze(F_trials(3,a,1)),'og','markersize',10);
        plot3(squeeze(F_trials(1,a,time_zero)),squeeze(F_trials(2,a,time_zero)),squeeze(F_trials(3,a,time_zero)),'.g','markersize',10);
        plot3(squeeze(F_trials(1,a,time_five)),squeeze(F_trials(2,a,time_five)),squeeze(F_trials(3,a,time_five)),'vg','markersize',10);
    end 
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'g','linewidth',2);
    plot3(mtrials(1,1),mtrials(2,1),mtrials(3,1),'og','markersize',10);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.g','markersize',10);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'vg','markersize',10);
    
    figure(456);
    subplot(3,1,1);
    shadedErrorBar(time_trial,mtrials(1,:),semtrials(1,:),'g');
    hold on
    subplot(3,1,2);
    shadedErrorBar(time_trial,mtrials(2,:),semtrials(2,:),'g');
    hold on
    subplot(3,1,3);
    shadedErrorBar(time_trial,mtrials(3,:),semtrials(3,:),'g');
    hold on
    
elseif strmatch('shock',label)
    figure(123);
    for a=1:length(trial_type)
        plot3(squeeze(F_trials(1,a,:)),squeeze(F_trials(2,a,:)),squeeze(F_trials(3,a,:)),'-r');
        hold on
        plot3(squeeze(F_trials(1,a,1)),squeeze(F_trials(2,a,1)),squeeze(F_trials(3,a,1)),'or','markersize',10);
        plot3(squeeze(F_trials(1,a,time_zero)),squeeze(F_trials(2,a,time_zero)),squeeze(F_trials(3,a,time_zero)),'.r','markersize',10);
        plot3(squeeze(F_trials(1,a,time_five)),squeeze(F_trials(2,a,time_five)),squeeze(F_trials(3,a,time_five)),'vr','markersize',10);
    end
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'r','linewidth',2);
    plot3(mtrials(1,1),mtrials(2,1),mtrials(3,1),'or','markersize',10);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.r','markersize',10);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'vr','markersize',10);
    
    figure(456);
    subplot(3,1,1);
    shadedErrorBar(time_trial,mtrials(1,:),semtrials(1,:),'r');
    hold on
    subplot(3,1,2);
    shadedErrorBar(time_trial,mtrials(2,:),semtrials(2,:),'r');
    hold on
    subplot(3,1,3);
    shadedErrorBar(time_trial,mtrials(3,:),semtrials(3,:),'r');
    hold on
    
    
elseif strmatch('miss',label)
    figure(123);
    for a=1:length(trial_type)
        plot3(squeeze(F_trials(1,a,:)),squeeze(F_trials(2,a,:)),squeeze(F_trials(3,a,:)),'-b');
        hold on
    end
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'b','linewidth',2);
    plot3(mtrials(1,1),mtrials(2,1),mtrials(3,1),'ob','markersize',20);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.b','markersize',20);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'vb','markersize',20);
    
    figure(456);
    subplot(3,1,1);
    shadedErrorBar(time_trial,mtrials(1,:),semtrials(1,:),'b');
    hold on
    subplot(3,1,2);
    shadedErrorBar(time_trial,mtrials(2,:),semtrials(2,:),'b');
    hold on
    subplot(3,1,3);
    shadedErrorBar(time_trial,mtrials(3,:),semtrials(3,:),'b');
    hold on
end

end