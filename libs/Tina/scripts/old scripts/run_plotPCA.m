function run_plotPCA(pathname,filedate,fn,frametime)

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
plot_pre=5;
plot_post=4;
pre_ind=round(plot_pre/frametime);
post_ind=round(plot_post/frametime);
total_samples=pre_ind+post_ind+1;
if plot_pre>0
    time_trial=time(1:total_samples)-time(pre_ind);
else
    time_trial=time(1:total_samples);
end

% create dfof for PCA
X=[];
for a=1:length(tone_start_2p)
    X=[X mean(dfof(:,tone_start_2p(a)-pre_ind:tone_start_2p(a)+post_ind),2)];
end

for a=1:size(X,1)
    neg=find(X(a,:)<0);
    X(a,neg)=0;
end

% % find active cells
% del=[];
% for a=1:size(X,1)
%     currcell=X(a,:);
%     threshold=6*mad(currcell,1);
%     active=[];
%     active=find(currcell>threshold);
%     if length(active)<10
%         del=[del a];
%     end
% end
% X(del,:)=[];

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

num_trials=length(tone_start_2p);

% plot trajectories aligned to the tone
get_traj(F,reward_trials,'reward');
get_traj(F,shock_trials,'shock');
get_traj(F,missed_trials,'miss');


if ~exist(strcat(currpath,'/results/'),'dir')
    mkdir(strcat(currpath,'/results/'));
end


savename=strcat(currpath,'/results/',fn,'_pca.mat');
save(savename,'time_trial','mpca_shock_tone','sempca_shock_tone','pca_trials_shock_tone','mpca_reward_tone','sempca_reward_tone','pca_trials_reward_tone','mpca_miss_tone','sempca_miss_tone','pca_trials_miss_tone');


end


function get_traj(F,trial_type,label)

dfof=F;

figure(123);
hold on
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

if strmatch('reward',label)
    for b=1:length(trial_type)
        a=trial_type(b);
        plot3(F(1,a),F(2,a),F(3,a),'.g','markersize',20);
        hold on
    end 
elseif strmatch('shock',label)
    for b=1:length(trial_type)
        a=trial_type(b);
        plot3(F(1,a),F(2,a),F(3,a),'.r','markersize',20);
        hold on
    end
elseif strmatch('miss',label)
    for b=1:length(trial_type)
        a=trial_type(b);
        plot3(F(1,a),F(2,a),F(3,a),'.b','markersize',20);
        hold on
    end
end

end