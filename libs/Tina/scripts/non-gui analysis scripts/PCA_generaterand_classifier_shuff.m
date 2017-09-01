function PCA_generaterand_classifier_shuff(pathname,filedate,fn)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'),'lever_retract','time_behav','time','reward_trials','shock_trials');
frametime=0.132423184;
dfof=dfof_sub;

% find behavioral stimulus starts
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

% randomly re-assign reward and shock trials
pressed=sort([shock_trials reward_trials]);
numshock=length(shock_trials);

% calculate all trajectories aligned to the press
[pca_trials_pressed]=get_traj(F,lever_retract_2p,pressed,time_trial,pre_ind,post_ind);

savename=strcat(currpath,'/results/',fn,'_pca_forshuff.mat');
save(savename,'time_trial','pca_trials_pressed','numshock');

end


function [F_trials]=get_traj(F,alignment,trial_type,time_trial,pre_ind,post_ind)

time_zero=find(time_trial>=1,1,'first');
time_five=find(time_trial>=5,1,'first');

F_trials=zeros(3,length(trial_type),length(time_trial));
for a=1:3
    for b=1:length(trial_type)
        curr_trial=trial_type(b);
        F_trials(a,b,:)=F(a,alignment(curr_trial)-pre_ind:alignment(curr_trial)+post_ind);
        F_trials(a,b,:)=smooth(F_trials(a,b,:),15);
    end
end

end