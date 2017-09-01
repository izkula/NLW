function run_plotLDA(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'));

dfof=dfof_sub;
time=[0:frametime:(length(dfof)-1)*frametime];

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

pressed_trials=sort([reward_trials shock_trials]);


% include this data in the PCA
plot_pre=0;
plot_post=5;
pre_ind=round(plot_pre/frametime);
post_ind=round(plot_post/frametime);
total_samples=pre_ind+post_ind+1;
if plot_pre>0
    time_trial=time(1:total_samples)-time(pre_ind);
else
    time_trial=time(1:total_samples);
end

% create dfof that only includes trials
X=[];
for a=1:length(tone_start_2p)
    X=[X dfof(:,tone_start_2p(a)-pre_ind:tone_start_2p(a)+post_ind)];
end
for a=1:size(X,1)
    maxval=max(X(a,:));
    X(a,:)=X(a,:)./maxval;
end

% find active cells
del=[];
for a=1:size(X,1)
    currcell=X(a,:);
    threshold=6*mad(currcell,1);
    
    currcell_trials=reshape(currcell,length(time_trial),length(tone_start_2p));
    active=[];
    for b=1:length(tone_start_2p)
        temp=[];
        currtrial=currcell_trials(:,b);
        temp=find(currtrial>threshold);
        if ~isempty(temp)
            active=[active b];
        end
    end
    
    if length(active)<.2*length(tone_start_2p)
        del=[del a];
    end
end
X(del,:)=[];

% create matrix for cov
X_all=zeros(length(tone_start_2p),length(time_trial),size(X,1));

for c=1:size(X,1)
    currcell=X(c,:);
    currcell_trials=reshape(currcell,length(time_trial),length(tone_start_2p));
    for a=1:length(tone_start_2p)
        curr_trial=currcell_trials(:,a)';
        X_all(a,:,c)=curr_trial;
    end
end

X_pressed=[];
for a=1:length(pressed_trials)
    curr_sig=squeeze(X_all(pressed_trials(a),:,:));
    X_pressed=[X_pressed; curr_sig];
end

X_missed=[];
for a=1:length(missed_trials)
    curr_sig=squeeze(X_all(missed_trials(a),:,:));
    X_missed=[X_missed; curr_sig];
end

cov_pressed=cov(X_pressed);
cov_missed=cov(X_missed);
mean_pressed=mean(X_pressed)';
mean_missed=mean(X_missed)';

w=(cov_pressed+cov_missed)\(mean_pressed-mean_missed);
w_cells=[w [1:length(w)]'];
w_cells=sortrows(w_cells);


pressed_traj=[];
for a=1:length(pressed_trials)
    current_sig=squeeze(X_all(pressed_trials(a),:,:));
    pressed_traj(a,:)=w'*current_sig';
end

missed_traj=[];
for a=1:length(missed_trials)
    current_sig=squeeze(X_all(missed_trials(a),:,:));
    missed_traj(a,:)=w'*current_sig';
end

figure(1);
for a=1:size(pressed_traj,1)
plot(time_trial,pressed_traj(a,:),'g');
hold on
end
for a=1:size(missed_traj,1)
plot(time_trial,missed_traj(a,:),'b');
hold on
end

figure(2);
shadedErrorBar(time_trial,mean(pressed_traj),std(pressed_traj)/sqrt(size(pressed_traj,1)),'g');
hold on
shadedErrorBar(time_trial,mean(missed_traj),std(missed_traj)/sqrt(size(missed_traj,1)),'b');
plot(time_trial,zeros(1,length(time_trial)),'-k');
xlim([0 5])



