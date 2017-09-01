close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

all_heatmap_pressed_pressed=[];
all_heatmap_pressed_missed=[];

all_heatmap_missed_pressed=[];
all_heatmap_missed_missed=[];

total_cells_pressed=[];
total_cells_missed=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/results/',filename,'.mat'));
    load(strcat(pathname,date,'/results/',filename,'_lasso_pressmiss.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_leveronly.mat'));
    
    t_0=find(time_trial>=0,1,'first');
    t_5=find(time_trial<=5,1,'last');
    
    num_cells=length(cells_lasso);
    missed_cells=find(weights_lasso>0);
    pressed_cells=find(weights_lasso<0);
 
    mtrials_missed_tone=mtrials_missed_tone(cells_lasso,:);
    mtrials_pressed_tone=mtrials_pressed_tone(cells_lasso,:);
    
    % normalize each cell by max value
    for a=1:size(mtrials_missed_tone,1)
        maxval_missed=max(mtrials_missed_tone(a,t_0:t_5));
        maxval_pressed=max(mtrials_pressed_tone(a,t_0:t_5));
        maxval=max([maxval_missed maxval_pressed]);
        mtrials_missed_tone(a,:)=mtrials_missed_tone(a,:)./maxval;
        mtrials_pressed_tone(a,:)=mtrials_pressed_tone(a,:)./maxval;
    end
        
    heatmap_pressed_pressed=mtrials_pressed_tone(pressed_cells,t_0:t_5);
    heatmap_pressed_missed=mtrials_missed_tone(pressed_cells,t_0:t_5);
    
    heatmap_missed_pressed=mtrials_pressed_tone(missed_cells,t_0:t_5);
    heatmap_missed_missed=mtrials_missed_tone(missed_cells,t_0:t_5);
    
    all_heatmap_pressed_pressed=[all_heatmap_pressed_pressed; heatmap_pressed_pressed];
    all_heatmap_pressed_missed=[all_heatmap_pressed_missed; heatmap_pressed_missed];
    
    all_heatmap_missed_pressed=[all_heatmap_missed_pressed; heatmap_missed_pressed];
    all_heatmap_missed_missed=[all_heatmap_missed_missed; heatmap_missed_missed];

    total_cells_pressed(z)=length(pressed_cells);
    total_cells_missed(z)=length(missed_cells);
    
    
%     missed_reward(z)=length(find(ismember(cells_lasso(missed_cells),reward_only)))/length(missed_cells);
%     missed_shock(z)=length(find(ismember(cells_lasso(missed_cells),shock_only)))/length(missed_cells);
%     pressed_reward(z)=length(find(ismember(cells_lasso(pressed_cells),reward_only)))/length(pressed_cells);
%     pressed_shock(z)=length(find(ismember(cells_lasso(pressed_cells),shock_only)))/length(pressed_cells);
    
%     missed_reward(z)=length(find(ismember(cells_lasso(missed_cells),reward_only)));
%     missed_shock(z)=length(find(ismember(cells_lasso(missed_cells),shock_only)));
%     pressed_reward(z)=length(find(ismember(cells_lasso(pressed_cells),reward_only)));
%     pressed_shock(z)=length(find(ismember(cells_lasso(pressed_cells),shock_only)));
end

border_pressed(1)=total_cells_pressed(1);
border_missed(1)=total_cells_missed(1);
for a=2:length(dates)
    border_pressed(a)=total_cells_pressed(a)+border_pressed(a-1);
    border_missed(a)=total_cells_missed(a)+border_missed(a-1);
end
border_pressed=border_pressed+0.5;
border_missed=border_missed+0.5;

figure;
subplot(1,2,1);
imagesc((all_heatmap_missed_missed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
subplot(1,2,2);
imagesc((all_heatmap_missed_pressed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);

figure;
subplot(1,2,1);
imagesc((all_heatmap_pressed_missed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_pressed' border_pressed'],'-w','linewidth',2);

subplot(1,2,2);
imagesc((all_heatmap_pressed_pressed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_pressed' border_pressed'],'-w','linewidth',2);

% figure;
% bar([1:4],[mean(missed_shock) mean(missed_reward) mean(pressed_shock) mean(pressed_reward)])
% hold on
% errorbar([1:4],[mean(missed_shock) mean(missed_reward) mean(pressed_shock) mean(pressed_reward)],[std(missed_shock) std(missed_reward) std(pressed_shock) std(pressed_reward)]./sqrt(5))

%%

all_heatmap_pressed_shock=[];
all_heatmap_pressed_reward=[];
all_heatmap_pressed_diff=[];

all_heatmap_missed_shock=[];
all_heatmap_missed_reward=[];
all_heatmap_missed_diff=[];

t_neg1=find(time_trial>=-1,1,'first');

for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/results/',filename,'.mat'));
    load(strcat(pathname,date,'/results/',filename,'_lasso_pressmiss.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_leveronly.mat'));
    
    t_1=find(time_trial>=1,1,'first');
    t_2=find(time_trial>=2,1,'first');
    t_3=find(time_trial>=3,1,'first');
    
    num_cells=length(cells_lasso);
    missed_cells=find(weights_lasso>0);
    pressed_cells=find(weights_lasso<0);

    mtrials_shock_press=mtrials_shock_press(cells_lasso,:);
    mtrials_reward_press=mtrials_reward_press(cells_lasso,:);
    
%     diff_shockreward=[];
%     for a=1:size(mtrials_shock_press,1)
%         meanval_shock=mean(mtrials_shock_press(a,t_1:t_2));
%         meanval_reward=mean(mtrials_reward_press(a,t_1:t_2));
%         diff_shockreward(a)=meanval_shock-meanval_reward;
%     end
%      
%     figure(541);
%     plot(weights_lasso,diff_shockreward,'o','color',colors(z,:));
%     hold on
  
%     % smooth data
%     for a=1:size(mtrials_shock_press,1)
%         mtrials_shock_press(a,:)=smooth(mtrials_shock_press(a,:),4);
%         mtrials_reward_press(a,:)=smooth(mtrials_reward_press(a,:),4);
%     end

%     % normalize each cell by max value
%     for a=1:size(mtrials_shock_press,1)
%         maxval=max(mtrials_shock_press(a,t_0:t_5));
%         mtrials_shock_press(a,:)=mtrials_shock_press(a,:)./maxval;
%     end
%     % normalize each cell by max value
%     for a=1:size(mtrials_reward_press,1)
%         maxval=max(mtrials_reward_press(a,t_0:t_5));
%         mtrials_reward_press(a,:)=mtrials_reward_press(a,:)./maxval;
%     end

    heatmap_pressed_shock=mtrials_shock_press(pressed_cells,t_neg1:t_5);
    heatmap_pressed_reward=mtrials_reward_press(pressed_cells,t_neg1:t_5);
    
%     % sort each cell by max response
%     indx=[];
%     for a=1:size(heatmap_pressed_shock,1)
%         [indx(a)]=mean(heatmap_pressed_shock(a,:));
%     end
%     order=[indx' [1:size(heatmap_pressed_shock,1)]'];
%     order=sortrows(order,-1);
%     order=order(:,2);
%     heatmap_pressed_shock=heatmap_pressed_shock(order,:);
%     indx=[];
%     for a=1:size(heatmap_pressed_reward,1)
%         [indx(a)]=mean(heatmap_pressed_reward(a,:));
%     end
%     order=[indx' [1:size(heatmap_pressed_reward,1)]'];
%     order=sortrows(order,-1);
%     order=order(:,2);
%     heatmap_pressed_reward=heatmap_pressed_reward(order,:);
%     
    heatmap_missed_shock=mtrials_shock_press(missed_cells,t_neg1:t_5);
    heatmap_missed_reward=mtrials_reward_press(missed_cells,t_neg1:t_5);
%     % sort each cell by max value
%     indx=[];
%     for a=1:size(heatmap_missed_shock,1)
%         [indx(a)]=mean(heatmap_missed_shock(a,:));
%     end
%     order=[indx' [1:size(heatmap_missed_shock,1)]'];
%     order=sortrows(order,-1);
%     order=order(:,2);
%     heatmap_missed_shock=heatmap_missed_shock(order,:);
%     indx=[];
%     for a=1:size(heatmap_missed_reward,1)
%         [indx(a)]=mean(heatmap_missed_reward(a,:));
%     end
%     order=[indx' [1:size(heatmap_missed_reward,1)]'];
%     order=sortrows(order,-1);
%     order=order(:,2);
%     heatmap_missed_reward=heatmap_missed_reward(order,:);
    
    
    
    all_heatmap_pressed_shock=[all_heatmap_pressed_shock; (heatmap_pressed_shock)];
    all_heatmap_pressed_reward=[all_heatmap_pressed_reward; (heatmap_pressed_reward)];
    all_heatmap_pressed_diff=[all_heatmap_pressed_diff; (heatmap_pressed_shock-heatmap_pressed_reward)];
    
    all_heatmap_missed_shock=[all_heatmap_missed_shock; (heatmap_missed_shock)];
    all_heatmap_missed_reward=[all_heatmap_missed_reward; (heatmap_missed_reward)];
    all_heatmap_missed_diff=[all_heatmap_missed_diff; (heatmap_missed_shock-heatmap_missed_reward)];
    
end

time_axis=time_trial(t_neg1:t_5);
t_lever=find(time_axis>=0,1,'first');
t_press=find(time_axis>=1,1,'first');
t_press2=find(time_axis>=2,1,'first');
t_press3=find(time_axis>=3,1,'first');

% smooth data
for a=1:size(all_heatmap_missed_shock,1)
    all_heatmap_missed_shock(a,:)=smooth(all_heatmap_missed_shock(a,:),4);
    all_heatmap_missed_reward(a,:)=smooth(all_heatmap_missed_reward(a,:),4);
end
for a=1:size(all_heatmap_pressed_shock,1)
    all_heatmap_pressed_shock(a,:)=smooth(all_heatmap_pressed_shock(a,:),4);
    all_heatmap_pressed_reward(a,:)=smooth(all_heatmap_pressed_reward(a,:),4);
end

% indx=[];
% for a=1:size(all_heatmap_missed_shock,1)
%     indx(a)=-mean(all_heatmap_missed_shock(a,1:t_press-1))+mean(all_heatmap_missed_shock(a,t_press:t_press3));
% end
% order=[indx' [1:size(all_heatmap_missed_shock,1)]'];
% order=sortrows(order,-1);
% order=order(:,2);
% all_heatmap_missed_shock=all_heatmap_missed_shock(order,:);
% 
% indx=[];
% for a=1:size(all_heatmap_missed_reward,1)
%     indx(a)=-mean(all_heatmap_missed_reward(a,1:t_press-1))+mean(all_heatmap_missed_reward(a,t_press:t_press3));
% end
% order=[indx' [1:size(all_heatmap_missed_reward,1)]'];
% order=sortrows(order,-1);
% order=order(:,2);
% all_heatmap_missed_reward=all_heatmap_missed_reward(order,:);
% 
% indx=[];
% for a=1:size(all_heatmap_pressed_reward,1)
%     indx(a)=-mean(all_heatmap_pressed_reward(a,1:t_press-1))+mean(all_heatmap_pressed_reward(a,t_press:t_press3));
% end
% order=[indx' [1:size(all_heatmap_pressed_reward,1)]'];
% order=sortrows(order,-1);
% order=order(:,2);
% all_heatmap_pressed_reward=all_heatmap_pressed_reward(order,:);
% 
% indx=[];
% for a=1:size(all_heatmap_pressed_shock,1)
%     indx(a)=-mean(all_heatmap_pressed_shock(a,1:t_press-1))+mean(all_heatmap_pressed_shock(a,t_press:t_press3));
% end
% order=[indx' [1:size(all_heatmap_pressed_shock,1)]'];
% order=sortrows(order,-1);
% order=order(:,2);
% all_heatmap_pressed_shock=all_heatmap_pressed_shock(order,:);

figure;
subplot(1,2,1);
imagesc((all_heatmap_missed_shock));
colorbar
caxis([0 1.5])
hold on
% plot([0.5 size(all_heatmap_pressed_shock,2)+0.5],[border_missed' border_missed'],'-w','linewidth',2);
%plot([t_lever t_lever],[0.5 size(all_heatmap_missed_shock,1)],'-w','linewidth',2);
plot([t_press t_press],[0.5 size(all_heatmap_missed_shock,1)],'-w','linewidth',2);
subplot(1,2,2);
imagesc((all_heatmap_missed_reward));
colorbar
caxis([0 1.5])
hold on
% plot([0.5 size(all_heatmap_pressed_shock,2)+0.5],[border_missed' border_missed'],'-w','linewidth',2);
% plot([t_lever t_lever],[0.5 size(all_heatmap_missed_shock,1)],'-w','linewidth',2);
plot([t_press t_press],[0.5 size(all_heatmap_missed_shock,1)],'-w','linewidth',2);
figure;
subplot(1,2,1);
imagesc((all_heatmap_pressed_shock));
colorbar
caxis([0 1.5])
hold on
% plot([0.5 size(all_heatmap_pressed_shock,2)+0.5],[border_pressed' border_pressed'],'-w','linewidth',2);
% plot([t_lever t_lever],[0.5 size(all_heatmap_pressed_shock,1)],'-w','linewidth',2);
plot([t_press t_press],[0.5 size(all_heatmap_pressed_shock,1)],'-w','linewidth',2);subplot(1,2,2);
imagesc((all_heatmap_pressed_reward));
colorbar
caxis([0 1.5])
hold on
% plot([0.5 size(all_heatmap_pressed_shock,2)+0.5],[border_pressed' border_pressed'],'-w','linewidth',2);
% plot([t_lever t_lever],[0.5 size(all_heatmap_pressed_shock,1)],'-w','linewidth',2);
plot([t_press t_press],[0.5 size(all_heatmap_pressed_shock,1)],'-w','linewidth',2);

figure;
shadedErrorBar(time_axis-1,mean(all_heatmap_missed_shock),std(all_heatmap_missed_shock)/sqrt(size(all_heatmap_missed_shock,1)),'r');
hold on
plot([0 0],[-0.1 0.6],'--k');
plot([-1 -1],[-0.1 0.6],'--k');
figure;
shadedErrorBar(time_axis-1,mean(all_heatmap_missed_reward),std(all_heatmap_missed_reward)/sqrt(size(all_heatmap_missed_reward,1)),'b');
hold on
plot([0 0],[-0.1 0.6],'--k');
plot([-1 -1],[-0.1 0.6],'--k');
figure;
shadedErrorBar(time_axis-1,mean(all_heatmap_pressed_shock),std(all_heatmap_pressed_shock)/sqrt(size(all_heatmap_pressed_shock,1)),'r');
hold on
plot([0 0],[-0.1 0.6],'--k');
plot([-1 -1],[-0.1 0.6],'--k');
figure;
shadedErrorBar(time_axis-1,mean(all_heatmap_pressed_reward),std(all_heatmap_pressed_reward)/sqrt(size(all_heatmap_pressed_reward,1)),'b');
hold on
plot([0 0],[-0.1 0.6],'--k');
plot([-1 -1],[-0.1 0.6],'--k');

% % all_heatmap_missed_diff=all_heatmap_missed_shock-all_heatmap_missed_reward;
% % figure;
% % shadedErrorBar(time_axis,mean(all_heatmap_missed_diff),std(all_heatmap_missed_diff)/sqrt(size(all_heatmap_missed_diff,1)),'k');
% 

t_0=find(time_axis>=0,1,'first');
t_press=find(time_axis>=1,1,'first');
t_press2=find(time_axis>=2,1,'first');
t_press3=find(time_axis>=3,1,'first');

startpoint=t_press;
endpoint=t_press2;

baseline_press_shock=mean(all_heatmap_pressed_shock(:,t_0:t_press-1),2);
response_press_shock=mean(all_heatmap_pressed_shock(:,startpoint:endpoint),2);
diff_press_shock=response_press_shock-baseline_press_shock;

baseline_press_reward=mean(all_heatmap_pressed_reward(:,t_0:t_press-1),2);
response_press_reward=mean(all_heatmap_pressed_reward(:,startpoint:endpoint),2);
diff_press_reward=response_press_reward-baseline_press_reward;

baseline_miss_shock=mean(all_heatmap_missed_shock(:,t_0:t_press-1),2);
response_miss_shock=mean(all_heatmap_missed_shock(:,startpoint:endpoint),2);
diff_miss_shock=response_miss_shock-baseline_miss_shock;

baseline_miss_reward=mean(all_heatmap_missed_reward(:,t_0:t_press-1),2);
response_miss_reward=mean(all_heatmap_missed_reward(:,startpoint:endpoint),2);
diff_miss_reward=response_miss_reward-baseline_miss_reward;

figure;
bar([1 2],[mean(diff_miss_shock) mean(diff_press_shock)])
hold on
errorbar([1 2],[mean(diff_miss_shock) mean(diff_press_shock)],[std(diff_miss_shock)/sqrt(length(diff_miss_shock)) std(diff_press_shock)/sqrt(length(diff_press_shock))])

figure;
bar([1 2],[mean(diff_miss_reward) mean(diff_press_reward)])
hold on
errorbar([1 2],[mean(diff_miss_reward) mean(diff_press_reward)],[std(diff_miss_reward)/sqrt(length(diff_miss_reward)) std(diff_press_reward)/sqrt(length(diff_press_reward))])