% close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% dates={'20160429','20160429','20160429','20160429','20160528'};
% filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

dates={'20160529','20160529','20160730','20160730','20160731'};
filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

all_heatmap_shock_missed=[];
all_heatmap_shock_pressed=[];

all_heatmap_reward_missed=[];
all_heatmap_reward_pressed=[];

total_cells_shock=[];
total_cells_reward=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/results/',filename,'.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_leveronly.mat'));
    
    t_0=find(time_trial>=0,1,'first');
    t_5=find(time_trial<5,1,'last');
    
%     del=sort([shock_only]);
%     shock_only=[1:size(mtrials_pressed_tone,1)];
%     shock_only(del)=[];

%     % smooth data
%     for a=1:size(mtrials_pressed_tone,1)
%         mtrials_pressed_tone(a,:)=smooth(mtrials_pressed_tone(a,:),4);
%         mtrials_missed_tone(a,:)=smooth(mtrials_missed_tone(a,:),4);
%     end

    heatmap_shock_pressed=mtrials_pressed_tone(shock_only,t_0:t_5);
    heatmap_reward_pressed=mtrials_pressed_tone(reward_only,t_0:t_5);
    
    heatmap_shock_missed=mtrials_missed_tone(shock_only,t_0:t_5);
    heatmap_reward_missed=mtrials_missed_tone(reward_only,t_0:t_5);
    
    all_heatmap_shock_missed=[all_heatmap_shock_missed; (heatmap_shock_missed)];
    all_heatmap_reward_missed=[all_heatmap_reward_missed; (heatmap_reward_missed)];
    
    all_heatmap_shock_pressed=[all_heatmap_shock_pressed; (heatmap_shock_pressed)];
    all_heatmap_reward_pressed=[all_heatmap_reward_pressed; (heatmap_reward_pressed)];
    
    total_cells_reward(z)=length(reward_only);
    total_cells_shock(z)=length(shock_only);
end


border_shock(1)=total_cells_shock(1);
border_reward(1)=total_cells_reward(1);
for a=2:length(dates)
    border_shock(a)=total_cells_shock(a)+border_shock(a-1);
    border_reward(a)=total_cells_reward(a)+border_reward(a-1);
end
border_shock=border_shock+0.5;
border_reward=border_reward+0.5;


% smooth data
for a=1:size(all_heatmap_shock_missed,1)
    all_heatmap_shock_missed(a,:)=smooth(all_heatmap_shock_missed(a,:),4);
    all_heatmap_shock_pressed(a,:)=smooth(all_heatmap_shock_pressed(a,:),4);
end
for a=1:size(all_heatmap_reward_missed,1)
    all_heatmap_reward_missed(a,:)=smooth(all_heatmap_reward_missed(a,:),4);
    all_heatmap_reward_pressed(a,:)=smooth(all_heatmap_reward_pressed(a,:),4);
end



time_axis=time_trial(t_0:t_5);

figure;
subplot(1,2,1);
imagesc((all_heatmap_shock_missed));
colorbar
hold on
plot([0.5 size(all_heatmap_shock_missed,2)+0.5],[border_shock' border_shock'],'-w','linewidth',2);
subplot(1,2,2);
imagesc((all_heatmap_shock_pressed));
colorbar
hold on
plot([0.5 size(all_heatmap_shock_pressed,2)+0.5],[border_shock' border_shock'],'-w','linewidth',2);

figure;
subplot(1,2,1);
imagesc((all_heatmap_reward_missed));
colorbar
hold on
plot([0.5 size(all_heatmap_reward_missed,2)+0.5],[border_reward' border_reward'],'-w','linewidth',2);
subplot(1,2,2);
imagesc((all_heatmap_reward_pressed));
colorbar
hold on
plot([0.5 size(all_heatmap_reward_pressed,2)+0.5],[border_reward' border_reward'],'-w','linewidth',2);

diff_shock=all_heatmap_shock_missed-all_heatmap_shock_pressed;
diff_reward=all_heatmap_reward_missed-all_heatmap_reward_pressed;
% diff_shock=mean(all_heatmap_shock_missed,2)/mean(all_heatmap_shock_pressed,2);
% diff_reward=mean(all_heatmap_reward_missed,2)/mean(all_heatmap_reward_pressed,2);

indx=[];
for a=1:size(diff_shock,1)
    indx(a)=mean(diff_shock(a,:));
end
order=[indx' [1:size(diff_shock,1)]'];
order=sortrows(order,-1);
order=order(:,2);
diff_shock=diff_shock(order,:);
boundary_shock=find(mean(diff_shock,2)<0,1,'first')-0.5;
boundary_shock1=find(mean(diff_shock,2)<0.1,1,'first')-0.5;
boundary_shock2=find(mean(diff_shock,2)<-0.1,1,'first')-0.5;

indx=[];
for a=1:size(diff_reward,1)
    indx(a)=mean(diff_reward(a,:));
end
order=[indx' [1:size(diff_reward,1)]'];
order=sortrows(order,-1);
order=order(:,2);
diff_reward=diff_reward(order,:);
boundary_reward=find(mean(diff_reward,2)<0,1,'first')-0.5;
boundary_reward1=find(mean(diff_reward,2)<0.1,1,'first')-0.5;
boundary_reward2=find(mean(diff_reward,2)<-0.1,1,'first')-0.5;

figure;
subplot(1,2,1);
imagesc((diff_shock));
colorbar
hold on
colormap(redblue);
%plot([0.5 size(diff_shock,2)+0.5],[boundary_shock boundary_shock],'-k');
plot([0.5 size(diff_shock,2)+0.5],[boundary_shock1 boundary_shock1],'-k');
plot([0.5 size(diff_shock,2)+0.5],[boundary_shock2 boundary_shock2],'-k');

subplot(1,2,2);
imagesc((diff_reward));
colorbar
hold on
colormap(redblue);
%plot([0.5 size(diff_reward,2)+0.5],[boundary_reward boundary_reward],'-k');
plot([0.5 size(diff_reward,2)+0.5],[boundary_reward1 boundary_reward1],'-k');
plot([0.5 size(diff_reward,2)+0.5],[boundary_reward2 boundary_reward2],'-k');


figure; bar([1 2],[mean(mean(all_heatmap_shock_missed,2)) mean(mean(all_heatmap_shock_pressed,2))])
hold on
errorbar([1 2],[mean(mean(all_heatmap_shock_missed,2)) mean(mean(all_heatmap_shock_pressed,2))],[std(mean(all_heatmap_shock_missed,2))/sqrt(size(all_heatmap_shock_missed,1)) std(mean(all_heatmap_shock_pressed,2))/sqrt(size(all_heatmap_shock_missed,1))])

figure; bar([1 2],[mean(mean(all_heatmap_reward_missed,2)) mean(mean(all_heatmap_reward_pressed,2))])
hold on
errorbar([1 2],[mean(mean(all_heatmap_reward_missed,2)) mean(mean(all_heatmap_reward_pressed,2))],[std(mean(all_heatmap_reward_missed,2))/sqrt(size(all_heatmap_reward_missed,1)) std(mean(all_heatmap_reward_pressed,2))/sqrt(size(all_heatmap_reward_missed,1))])
