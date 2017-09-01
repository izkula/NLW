function run_plotlevertrials(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/behavior/',fn,'.mat'));

dfof=dfof_sub;
%time=[0:frametime:(length(dfof)-1)*frametime];
    
for a=1:length(tone_start)
    tone_start_2p(a)=find(time<=time_behav(tone_start(a)),1,'last');
end

for a=1:length(lever_extend)
    lever_extend_2p(a)=find(time<=time_behav(lever_extend(a)),1,'last');
end

for a=1:length(lever_retract)
    lever_retract_2p(a)=find(time<=time_behav(lever_retract(a)),1,'last');
end

plot_pre=2;
plot_post=5;
pre_ind=round(plot_pre/frametime);
post_ind=round(plot_post/frametime);
total_samples=pre_ind+post_ind+1;
time_trial=time(1:total_samples)-time(pre_ind);

pressed_trials=[reward_trials shock_trials];
pressed_trials=sort(pressed_trials);

% plot data aligned to the lever press
[mtrials_shock_press, semtrials_shock_press, dfof_trials_shock_press,order_shock,maxval_shock]=get_tone(dfof,time_trial,shock_trials,lever_retract_2p,pre_ind,post_ind,[],[],'shock');
[mtrials_reward_press, semtrials_reward_press, dfof_trials_reward_press,order_reward,maxval_reward]=get_tone(dfof,time_trial,reward_trials,lever_retract_2p,pre_ind,post_ind,order_shock,[],'reward');
%[mtrials_shock_press, semtrials_shock_press, dfof_trials_shock_press,order_shock,maxval_shock]=get_tone(dfof,time_trial,shock_trials,lever_retract_2p,pre_ind,post_ind,[],[],'shock');


% plot data aligned to the tone
[mtrials_pressed_tone, semtrials_pressed_tone, dfof_trials_pressed_tone,order_pressed,maxval_pressed]=get_tone(dfof,time_trial,pressed_trials,tone_start_2p,pre_ind,post_ind,[],[],'tone');
[mtrials_missed_tone, semtrials_missed_tone, dfof_trials_missed_tone,order_missed,maxval_missed]=get_tone(dfof,time_trial,missed_trials,tone_start_2p,pre_ind,post_ind,order_pressed,[],'missed');
% %[mtrials_missed_tone, semtrials_missed_tone, dfof_trials_missed_tone,order_missed,maxval_missed]=get_tone(dfof,time_trial,missed_trials,tone_start_2p,pre_ind,post_ind,[],maxval_pressed,'missed');

if ~exist(strcat(currpath,'/results/'),'dir')
    mkdir(strcat(currpath,'/results/'));
end


savename=strcat(currpath,'/results/',fn,'.mat');
save(savename,'time_trial','mtrials_shock_press','semtrials_shock_press','dfof_trials_shock_press','mtrials_reward_press','semtrials_reward_press','dfof_trials_reward_press','mtrials_pressed_tone','semtrials_pressed_tone','dfof_trials_pressed_tone','mtrials_missed_tone','semtrials_missed_tone','dfof_trials_missed_tone');


end


function [mtrials,semtrials,dfof_trials,order,maxval]=get_tone(dfof,time_trial,trial_type,tone,pre,post,order1,maxval1,label)

%figure;
dfof_trials=zeros(size(dfof,1),length(trial_type),length(time_trial));
for a=1:size(dfof,1)
    for b=1:length(trial_type)
        curr_trial=trial_type(b);
        dfof_trials(a,b,:)=dfof(a,tone(curr_trial)-pre:tone(curr_trial)+post);
    end
    mtrials(a,:)=mean(squeeze(dfof_trials(a,:,:)));
    semtrials(a,:)=std(squeeze(dfof_trials(a,:,:)))/sqrt(size(dfof_trials,2));
    
    %shadedErrorBar(time_trial,mtrials(a,:)+a*0.5,semtrials(a,:));
    %hold on
end

%plot([0 0],[0 size(dfof,1)/2],'k','linewidth',2);
%plot([5 5],[0 size(dfof,1)/2],'k','linewidth',2);


if ~isempty(order1)
    order=order1;
else 
    % determine order
    order=zeros(size(dfof,1),2);
    order(:,1)=[1:size(dfof,1)];
    for a=1:size(dfof,1)
        [maxval,position]=max(mtrials(a,:));
        order(a,2)=position;
    end

    order=sortrows(order,2);
    peakpos=order(:,2);
    order=order(:,1);
end

heatmap=zeros(length(order),length(time_trial));

if ~isempty(maxval1)
    maxval=maxval1;
    for a=1:length(order)
        %heatmap(a,:)=mtrials(order(a),:)./maxval(a);
        heatmap(a,:)=mtrials(order(a),:);
        %heatmap(a,:)=smooth(heatmap(a,:),3);
    end
else
    for a=1:length(order)
        maxval(a)=max(mtrials(order(a),:));
        %heatmap(a,:)=mtrials(order(a),:)./maxval(a);
        heatmap(a,:)=mtrials(order(a),:);
        %heatmap(a,:)=smooth(heatmap(a,:),3);
    end
end

time_zero=find(time_trial>=0,1,'first');
time_five=find(time_trial>=5,1,'first');
time_one=find(time_trial>=1,1,'first');

figure;
imagesc(heatmap);
colormap(jet);
% caxis([0 2]);
colorbar;
hold on

if strmatch('reward',label)
    plot([time_zero time_zero],[0 size(dfof,1)],'k','linewidth',2);
    plot([time_one time_one],[0 size(dfof,1)],'k','linewidth',2);
    title('Reward');
elseif strmatch('tone',label)
    plot([time_zero time_zero],[0 size(dfof,1)],'k','linewidth',2);
    plot([time_five time_five],[0 size(dfof,1)],'k','linewidth',2);
    title('Tone');
elseif strmatch('shock',label)
    plot([time_zero time_zero],[0 size(dfof,1)],'k','linewidth',2);
    plot([time_one time_one],[0 size(dfof,1)],'k','linewidth',2);
    title('Shock');
elseif strmatch('missed',label)
    plot([time_zero time_zero],[0 size(dfof,1)],'k','linewidth',2);
    plot([time_five time_five],[0 size(dfof,1)],'k','linewidth',2);
    title('Missed');
end

end


% % determine order
% order=zeros(size(dfof_filt,1),2);
% order(:,1)=[1:size(dfof_filt,1)];
% for a=1:size(dfof_filt,1)
%     [maxval,position]=max(mtrials(a,:));
%     order(a,2)=position;
% end
% 
% order=sortrows(order,2);
% %maxval=max(max(mtrials));
% 
% heatmap=zeros(length(order),length(time_trial));
% for a=1:length(order)
%     maxval=max(mtrials(order(a),:));
%     heatmap(a,:)=mtrials(order(a),:)./maxval;
%     heatmap(a,:)=smooth(heatmap(a,:),3);
% end
%     
% figure;
% imagesc(heatmap);
% %colormap(jet);
% caxis([0 1]);
% hold on
% plot([time_zero time_zero],[0 size(dfof_filt,1)],'k','linewidth',2);
% plot([time_five time_five],[0 size(dfof_filt,1)],'k','linewidth',2);
% 
% figure(500);
% for a=1:size(dfof_filt,1)
%    plot(time,dfof_filt(a,:)+(a-1));
%    hold on
% end
% lever_times=ones(1,length(time))*-0.5;
% lever_times(lever_start_2p)=size(dfof_filt,1);
% lever_times_end=ones(1,length(time))*-0.5;
% lever_times_end(lever_start_2p+find(time>=5,1,'first'))=size(dfof_filt,1);
% hold on
% plot(time,lever_times);
% plot(time,lever_times_end);
% 
% 
% %timelabels=time_trial([1 5 10 15 20 25 30 35 40]);
% %set(gca,'XTickLabel',timelabels)
% %
% % %%
% % 
% % alltrials=[];
% % for a=1:size(dfof_filt,1)
% %     trials=[];
% %     for b=1:length(lever_start)
% %         trials=[trials dfof_filt(a,lever_start_2p(b)-pre:lever_start_2p(b)+post)];
% %     end
% %     alltrials=[alltrials; trials];
% % end
% % 
% % % trialstarts=[1:40:40*31];
% % % trialends=trialstarts+39;
% % 
% % tina=zeros(size(dfof_filt,1),length(dfof_filt));
% % for a=1:size(dfof_filt,1)
% %     tina(a,:)=smooth(dfof_filt(a,:),10);
% % end
% % 
% % 
% % [F,pcs,eigs,var,pvar]=pca(tina);
% % 
% % figure(500);
% % plot(pcs(:,1));
% % hold on
% % plot(pcs(:,2),'r');
% % plot(pcs(:,3),'g');
% % legend('PC1','PC2','PC3');
% % figure(501);
% % plot(pvar);
% % 
% % figure;
% % for a=1:length(lever_start_2p)
% %     axis1=F(lever_start_2p(a)-pre:lever_start_2p(a)+post,1);
% %     axis2=F(lever_start_2p(a)-pre:lever_start_2p(a)+post,2);
% %     axis3=F(lever_start_2p(a)-pre:lever_start_2p(a)+post,3);
% %     plot3(smooth(axis1,3),smooth(axis2,3),smooth(axis3,3));
% %     hold on
% %     plot3(axis1(1),axis2(1),axis3(1),'.','markersize',30);
% %     %plot3(axis1(end),axis2(end),axis3(end),'>','markersize',10);
% %     %plot3(axis1(11),axis2(11),axis3(11),'.','markersize',30);
% %     %plot3(axis1(17),axis2(17),axis3(17),'^','markersize',10);
% % end
% % 
% % PC1=[];
% % PC2=[];
% % PC3=[];
% % for a=1:length(lever_start_2p)
% %     PC1=[PC1; F(lever_start_2p(a)-pre:lever_start_2p(a)+post,1)'];
% %     PC2=[PC2; F(lever_start_2p(a)-pre:lever_start_2p(a)+post,2)'];
% %     PC3=[PC3; F(lever_start_2p(a)-pre:lever_start_2p(a)+post,3)'];
% % end
% % 
% % distances=zeros(size(PC1,1),size(PC1,1),40);
% % for a=1:40
% %     distances(:,:,a)=squareform(pdist([PC1(:,a) PC2(:,a) PC3(:,a)]));
% %     
% % end
% % alldist=mean(distances,3);
% % 
% % numclusters=20;
% % IDX=kmedoids(alldist,numclusters);
% % 
% % for a=1:numclusters
% %     figure(a+200);
% %     clustertrials=find(IDX==a);
% %     for b=1:length(clustertrials)
% %         axis1=F(lever_start_2p(clustertrials(b))-pre:lever_start_2p(clustertrials(b))+post,1);
% %         axis2=F(lever_start_2p(clustertrials(b))-pre:lever_start_2p(clustertrials(b))+post,2);
% %         axis3=F(lever_start_2p(clustertrials(b))-pre:lever_start_2p(clustertrials(b))+post,3);
% %         plot3(smooth(axis1,3),smooth(axis2,3),smooth(axis3,3));
% %         hold on
% %         plot3(axis1(1),axis2(1),axis3(1),'.','markersize',20);
% %         plot3(axis1(end),axis2(end),axis3(end),'+','markersize',20);
% %     end
% % end
% % 
% % % try sorting by lick time?
% % 
% % % try sorting by the time elapsed since the previous lever press?
% % 
% % 
% % 
% % 
% % 
% % 
% % 