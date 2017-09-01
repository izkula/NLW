function run_plotLDA(pathname,filedate,fn,frametime)

load(strcat(pathname,filedate,'/results/',filename,'.mat'));

num_cells=size(dfof_trials_pressed_tone,1);
num_pressed=size(dfof_trials_pressed_tone,2);
num_missed=size(dfof_trials_missed_tone,2);
num_trials=num_pressed+num_missed;

t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial<=5,1,'last');

% create "X" variable for regression, only for t=0 to 5 of tone.
pressed=dfof_trials_pressed_tone(:,:,t_0:t_5);
pressed=reshape(pressed,[num_cells,size(pressed,2)*size(pressed,3)]);
missed=dfof_trials_missed_tone(:,:,t_0:t_5);
missed=reshape(missed,[num_cells,size(missed,2)*size(missed,3)]);
X=[pressed'; missed'];

% create matrix for cov
X_all=zeros(length(lever_retract_2p),length([t_0:t_5]),size(X,1));

for c=1:size(X,1)
    currcell=X(c,:);
    currcell_trials=reshape(currcell,length([t_0:t_5]),num_trials);
    for a=1:num_trials
        curr_trial=currcell_trials(:,a)';
        X_all(a,:,c)=curr_trial;
    end
end

% perform jackknife cross-validation
pressed_traj_test=[];
missed_traj_test=[];
w_cells_all=[];
for z=1:num_trials
    test_trial=z;
    train_missed=[];
    train_pressed=[];
    
    if ismember(test_trial,pressed_trials)
        test_type='pressed';
        del=find(pressed_trials==test_trial);
        train_missed=[1:num_missed];
        train_pressed=[1:num_pressed];
        train_pressed(del)=[];
    else
        test_type='missed';
        del=find(missed_trials==test_trial);
        train_missed=[1:num_missed];
        train_pressed=[1:num_pressed];
        train_missed(del)=[];
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
%     w_cells(:,1)=w_cells(:,1)/abs(max(w_cells(:,1)));
    w_cells_all(:,:,z)=w_cells;

%     for a=1:length(train_pressed)
%         curr_trial=pressed_trials(train_pressed(a));
%         current_sig=squeeze(X_all(curr_trial,:,:));
%         pressed_traj_train=[pressed_traj_train; w'*current_sig'];
%     end
% 
%     for a=1:length(train_missed)
%         curr_trial=missed_trials(train_missed(a));
%         current_sig=squeeze(X_all(curr_trial,:,:));
%         missed_traj_train=[missed_traj_train; w'*current_sig'];
%     end
    
    if strmatch('pressed',test_type)
        current_sig=squeeze(X_all(test_trial,:,:));
        pressed_traj_test=[pressed_traj_test; w'*current_sig'];
    else
        current_sig=squeeze(X_all(test_trial,:,:));
        missed_traj_test=[missed_traj_test; w'*current_sig'];
    end
end

figure(1);
%subplot(2,1,1);
for a=1:num_pressed
    %plot(time_trial,pressed_traj_test(a,:),'--g');
    %hold on
    pressed_traj_test(a,:)=smooth(pressed_traj_test(a,:),10);
    pressed_test(a)=mean(pressed_traj_test(a,:));
    %search=sort(pressed_traj_test(a,:),2,'descend');
%     pressed_test(a)=mean(pressed_traj_test(a,1:20));
end
shadedErrorBar(time_trial,mean(pressed_traj_test),std(pressed_traj_test)/sqrt(size(pressed_traj_test,1)),'g');
hold on
for a=1:num_missed
    %plot(time_trial,missed_traj_test(a,:),'--b');
    %hold on
    missed_traj_test(a,:)=smooth(missed_traj_test(a,:),10);
    missed_test(a)=mean(missed_traj_test(a,:));
    %search=sort(missed_traj_test(a,:),2,'ascend');
    %missed_test(a)=mean(missed_traj_test(a,1:20));
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

[Xauc,Yauc,~,AUC]=perfcurve(trial_types,[pressed_test missed_test],'pressed');

figure(2);
plot(Xauc,Yauc);
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

p_auc=[];

AUC_shuff=zeros(1,100);
for s=1:1000
    pressed_traj_test=[];
    missed_traj_test=[];
%     pressed_traj_train=[];
%     missed_traj_train=[];
    
    order=randperm(num_trials);
    pressed_trials_shuff=order(1:num_pressed);
    pressed_trials_shuff=sort(pressed_trials_shuff);
    missed_trials_shuff=order(num_pressed+1:end);
    missed_trials_shuff=sort(missed_trials_shuff);
    
    for z=1:length(tone_start_2p)
        
        test_trial=z;
        train_missed=[];
        train_pressed=[];
        
        if ismember(test_trial,pressed_trials_shuff)
            test_type='pressed';
            del=find(pressed_trials_shuff==test_trial);
            train_missed=[1:num_missed];
            train_pressed=[1:num_pressed];
            train_pressed(del)=[];
        else
            test_type='missed';
            del=find(missed_trials_shuff==test_trial);
            train_missed=[1:num_missed];
            train_pressed=[1:num_pressed];
            train_missed(del)=[];
        end
        
        X_missed_train=[];
        for a=1:length(train_missed)
            curr_trial=missed_trials_shuff(train_missed(a));
            curr_sig=squeeze(X_all(curr_trial,:,:));
            X_missed_train=[X_missed_train; curr_sig];
        end
        
        X_pressed_train=[];
        for a=1:length(train_pressed)
            curr_trial=pressed_trials_shuff(train_pressed(a));
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
         
%         for a=1:length(train_pressed)
%             curr_trial=pressed_trials_shuff(train_pressed(a));
%             current_sig=squeeze(X_all(curr_trial,:,:));
%             pressed_traj_train=[pressed_traj_train; w'*current_sig'];
%         end
%         
%         for a=1:length(train_missed)
%             curr_trial=missed_trials_shuff(train_missed(a));
%             current_sig=squeeze(X_all(curr_trial,:,:));
%             missed_traj_train=[missed_traj_train; w'*current_sig'];
%         end
        
        if strmatch('pressed',test_type)
            current_sig=squeeze(X_all(test_trial,:,:));
            pressed_traj_test=[pressed_traj_test; w'*current_sig'];
        else
            current_sig=squeeze(X_all(test_trial,:,:));
            missed_traj_test=[missed_traj_test; w'*current_sig'];
        end
    end
    
    pressed_test=[];
    for a=1:num_pressed
        pressed_traj_test(a,:)=smooth(pressed_traj_test(a,:),10);
        pressed_test(a)=mean(pressed_traj_test(a,:));
        %search=sort(pressed_traj_test(a,:),2,'descend');
%         pressed_test(a)=mean(pressed_traj_test(a,1:20));
    end
    missed_test=[];
    for a=1:num_missed
        missed_traj_test(a,:)=smooth(missed_traj_test(a,:),10);
        missed_test(a)=mean(missed_traj_test(a,:));
        %search=sort(missed_traj_test(a,:),2,'ascend');
        %         missed_test(a)=mean(missed_traj_test(a,1:20));
    end
    
    trial_types_shuff=[];
    for a=1:num_pressed+num_missed
        if a<=num_pressed
            trial_types_shuff{a}='pressed';
        else
            trial_types_shuff{a}='missed';
        end
    end
    
    [~,~,~,AUC_shuff(s)]=perfcurve(trial_types_shuff,[pressed_test missed_test],'pressed');
    display(s);
end

fake=length(find(AUC_shuff>=AUC));
p_auc=fake/100;

display(numcells);
display(p_auc);
display(AUC);

save(savename,'AUC','top5_missed','top5_pressed','p_auc');