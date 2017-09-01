function run_plotLDA_5fold(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'));

dfof=dfof_sub;
%time=[0:frametime:(length(dfof)-1)*frametime];

%flag='tone';
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

% include this data in the LDA
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

if strmatch('tone',flag,'exact')
    pressed_trials=sort([reward_trials shock_trials]);
    X=[];
    for a=1:length(tone_start_2p)
        X=[X dfof(:,tone_start_2p(a)-pre_ind:tone_start_2p(a)+post_ind)];
    end
    savename=strcat(currpath,'/results/',fn,'_lda_tone.mat');
else
    pressed_trials=sort([reward_trials shock_trials]);
    X=[];
    for a=1:length(lever_retract_2p)
        X=[X dfof(:,lever_retract_2p(a)-pre_ind:lever_retract_2p(a)+post_ind)];
    end
    savename=strcat(currpath,'/results/',fn,'_lda_press.mat');
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
display(size(X,1));

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

% perform 5-fold cross-validation
num_missed=length(missed_trials);
num_pressed=length(pressed_trials);

bin_missed=[1:floor(num_missed/5):num_missed];
bin_missed=bin_missed(1:5);
binsize_missed=floor(num_missed/5);
bin_pressed=[1:floor(num_pressed/5):num_pressed];
bin_pressed=bin_pressed(1:5);
binsize_pressed=floor(num_pressed/5);


pressed_traj_test=[];
missed_traj_test=[];
pressed_traj_train=[];
missed_traj_train=[];
w_cells_all=[];
for z=1:5
    
    if z==5
        test_missed=[bin_missed(z):num_missed];
        train_missed=[1:num_missed];
        train_missed(test_missed)=[];

        test_pressed=[bin_pressed(z):num_pressed];
        train_pressed=[1:num_pressed];
        train_pressed(test_pressed)=[];
    else
        test_missed=[bin_missed(z):bin_missed(z)+binsize_missed-1];
        train_missed=[1:num_missed];
        train_missed(test_missed)=[];

        test_pressed=[bin_pressed(z):bin_pressed(z)+binsize_pressed-1];
        train_pressed=[1:num_pressed];
        train_pressed(test_pressed)=[];
    end
    
    X_missed_train=[];
    for a=1:length(train_missed)
        curr_trial=missed_trials(train_missed(a));
        curr_sig=squeeze(X_all(curr_trial,:,:));
        X_missed_train=[X_missed_train; curr_sig];
    end
    
    X_pressed_train=[];
    for a=1:length(train_pressed)
        curr_trial=pressed_trials(train_pressed(a));
        curr_sig=squeeze(X_all(curr_trial,:,:));
        X_pressed_train=[X_pressed_train; curr_sig];
    end

    cov_pressed=cov(X_pressed_train);
    cov_missed=cov(X_missed_train);
    mean_pressed=mean(X_pressed_train)';
    mean_missed=mean(X_missed_train)';

    w=(cov_pressed+cov_missed)\(mean_pressed-mean_missed);
    w_cells=[w [1:length(w)]'];
    w_cells=sortrows(w_cells);
    w_cells_all(:,:,z)=w_cells;
    
    for a=1:length(train_pressed)
        curr_trial=pressed_trials(train_pressed(a));
        current_sig=squeeze(X_all(curr_trial,:,:));
        pressed_traj_train=[pressed_traj_train; w'*current_sig'];
    end

    for a=1:length(train_missed)
        curr_trial=missed_trials(train_missed(a));
        current_sig=squeeze(X_all(curr_trial,:,:));
        missed_traj_train=[missed_traj_train; w'*current_sig'];
    end

    for a=1:length(test_pressed)
        curr_trial=pressed_trials(test_pressed(a));
        current_sig=squeeze(X_all(curr_trial,:,:));
        pressed_traj_test=[pressed_traj_test; w'*current_sig'];
    end

    for a=1:length(test_missed)
        curr_trial=missed_trials(test_missed(a));
        current_sig=squeeze(X_all(curr_trial,:,:));
        missed_traj_test=[missed_traj_test; w'*current_sig'];
    end
    
end

figure(3);
%subplot(2,1,1);
for a=1:num_pressed
    %plot(time_trial,pressed_traj_test(a,:),'--g');
    %hold on
    pressed_traj_test(a,:)=smooth(pressed_traj_test(a,:),10);
    pressed_test(a)=mean(pressed_traj_test(a,:));
%     pressed_test(a)=max(pressed_traj_test(a,:));
end
shadedErrorBar(time_trial,mean(pressed_traj_test),std(pressed_traj_test)/sqrt(size(pressed_traj_test,1)),'g');
hold on
for a=1:num_missed
    %plot(time_trial,missed_traj_test(a,:),'--b');
    %hold on
    missed_traj_test(a,:)=smooth(missed_traj_test(a,:),10);
    missed_test(a)=mean(missed_traj_test(a,:));
%     missed_test(a)=min(missed_traj_test(a,:));
end
shadedErrorBar(time_trial,mean(missed_traj_test),std(missed_traj_test)/sqrt(size(missed_traj_test,1)),'b');
plot(time_trial,zeros(1,length(time_trial)),'-k');
hold off
% subplot(2,1,2);
% plot([1:num_pressed],pressed_test,'.-g');
% hold on
% plot([1:num_missed],missed_test,'.-b');
% plot([1:max([num_pressed num_missed])],zeros(1,max([num_pressed num_missed])),'-k');
% hold off

for a=1:num_pressed+num_missed
    if a<=num_pressed
        trial_types{a}='pressed';
    else
        trial_types{a}='missed';
    end
end

[X,Y,T,AUC]=perfcurve(trial_types,[pressed_test missed_test],'pressed');

figure(4);
plot(X,Y);
xlabel('False positive rate');
ylabel('True positive rate');
text(0.9,0.1,num2str(AUC),'HorizontalAlignment','c');
hold off

missed_cells=w_cells_all(end-4:end,2,:);
missed_cells=reshape(missed_cells,1,size(missed_cells,1)*size(missed_cells,3));
figure;
missed_hist=histogram(missed_cells,'binwidth',1);
missed_counts=[missed_hist.Values' [1:length(missed_hist.Values)]'];
missed_counts=sortrows(missed_counts,1);
top5_missed_inds=missed_counts(end-4:end,2);
top5_missed=missed_hist.BinEdges(top5_missed_inds);

pressed_cells=w_cells_all(1:5,2,:);
pressed_cells=reshape(pressed_cells,1,size(pressed_cells,1)*size(pressed_cells,3));
figure;
pressed_hist=histogram(pressed_cells,'binwidth',1);
pressed_counts=[pressed_hist.Values' [1:length(pressed_hist.Values)]'];
pressed_counts=sortrows(pressed_counts,1);
top5_pressed_inds=pressed_counts(end-4:end,2);
top5_pressed=pressed_hist.BinEdges(top5_pressed_inds);

save(savename,'AUC','top5_missed','top5_pressed');