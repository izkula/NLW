function run_plotPCA_old(pathname,filedate,fn,frametime)

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

% include this data in the PCA
plot_pre=2;
plot_post=10;
pre_ind=round(plot_pre/frametime);
post_ind=round(plot_post/frametime);
total_samples=pre_ind+post_ind+1;
time_trial=time(1:total_samples)-time(pre_ind);

pressed_trials=[reward_trials shock_trials];
pressed_trials=sort(pressed_trials);

% create dfof for PCA
X=[];
for a=1:length(tone_start_2p(pressed_trials))
    X=[X dfof(:,tone_start_2p(pressed_trials(a))-pre_ind:tone_start_2p(pressed_trials(a))+post_ind)];
end
for a=1:size(X,1)
    maxval=max(X(a,:));
    X(a,:)=X(a,:)./maxval;
end

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

num_trials=length(tone_start_2p(pressed_trials));

% plot trajectories aligned to the tone
[mpca_reward_tone, sempca_reward_tone, pca_trials_reward_tone]=get_traj(F,num_trials,time_trial,pressed_trials(ismember(pressed_trials,reward_trials)),'reward');
[mpca_shock_tone, sempca_shock_tone, pca_trials_shock_tone]=get_traj(F,num_trials,time_trial,pressed_trials(ismember(pressed_trials,shock_trials)),'shock');
%[mpca_miss_tone, sempca_miss_tone, pca_trials_miss_tone]=get_traj(F,num_trials,time_trial,missed_trials,'miss');


if ~exist(strcat(currpath,'/results/'),'dir')
    mkdir(strcat(currpath,'/results/'));
end


savename=strcat(currpath,'/results/',fn,'_pca.mat');
%save(savename,'time_trial','mpca_shock_tone','sempca_shock_tone','pca_trials_shock_tone','mpca_reward_tone','sempca_reward_tone','pca_trials_reward_tone','mpca_miss_tone','sempca_miss_tone','pca_trials_miss_tone');
save(savename,'time_trial','mpca_shock_tone','sempca_shock_tone','pca_trials_shock_tone','mpca_reward_tone','sempca_reward_tone','pca_trials_reward_tone');


end


function [mtrials,semtrials,dfof_trials]=get_traj(F,num_trials,time_trial,trial_type,label)

axis1=F(1,:);
axis1=reshape(axis1,num_trials,length(time_trial));
axis2=F(2,:);
axis2=reshape(axis2,num_trials,length(time_trial));
axis3=F(3,:);
axis3=reshape(axis3,num_trials,length(time_trial));

dfof(1,:,:)=axis1;
dfof(2,:,:)=axis2;
dfof(3,:,:)=axis3;

time_zero=find(time_trial>=0,1,'first');
time_five=find(time_trial>=5,1,'first');

dfof_trials=zeros(3,length(trial_type),length(time_trial));
for a=1:3
    for b=1:length(trial_type)
        curr_trial=trial_type(b);
        dfof_trials(a,b,:)=dfof(a,b,:);
        dfof_trials(a,b,:)=smooth(dfof_trials(a,b,:),15);
    end
    mtrials(a,:)=mean(squeeze(dfof_trials(a,:,:)));
    semtrials(a,:)=std(squeeze(dfof_trials(a,:,:)))/sqrt(size(dfof_trials,2));
end

figure(123);

if strmatch('reward',label)
    for a=1:length(trial_type)
        %plot3(squeeze(dfof_trials(1,a,:)),squeeze(dfof_trials(2,a,:)),squeeze(dfof_trials(3,a,:)),'g');
        hold on
    end 
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'g','linewidth',2);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.g','markersize',30);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'xg','markersize',20);
elseif strmatch('shock',label)
    for a=1:length(trial_type)
        %plot3(squeeze(dfof_trials(1,a,:)),squeeze(dfof_trials(2,a,:)),squeeze(dfof_trials(3,a,:)),'r');
        hold on
    end
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'r','linewidth',2);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.r','markersize',30);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'xr','markersize',20);
elseif strmatch('miss',label)
    for a=1:length(trial_type)
        %plot3(squeeze(dfof_trials(1,a,:)),squeeze(dfof_trials(2,a,:)),squeeze(dfof_trials(3,a,:)),'b');
        hold on
    end
    plot3(mtrials(1,:),mtrials(2,:),mtrials(3,:),'b','linewidth',2);
    plot3(mtrials(1,time_zero),mtrials(2,time_zero),mtrials(3,time_zero),'.b','markersize',30);
    plot3(mtrials(1,time_five),mtrials(2,time_five),mtrials(3,time_five),'xb','markersize',20);
end

end