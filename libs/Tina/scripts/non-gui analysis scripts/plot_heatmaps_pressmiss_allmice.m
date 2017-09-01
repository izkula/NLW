clear all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';


dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

all_mtrials_pressed=[];
all_mtrials_missed=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/results/',filename,'.mat'));
        
    temp=squeeze(mean(dfof_trials_pressed_tone,2));
    all_mtrials_pressed=[all_mtrials_pressed; temp(shock_only,:)];
    temp=squeeze(mean(dfof_trials_missed_tone,2));
    all_mtrials_missed=[all_mtrials_missed; temp(shock_only,:)];

end

t_0=find(time_trial>=0,1,'first');
t_1=find(time_trial>=1,1,'first');
t_2=find(time_trial>=2,1,'first');
t_3=find(time_trial>=3,1,'first');
t_4=find(time_trial>=4,1,'first');
t_5=find(time_trial<5,1,'last');

num_cells=size(all_mtrials_missed,1);

ind=[];
for a=1:num_cells
    all_mtrials_missed(a,:)=smooth(all_mtrials_missed(a,:),4);
    all_mtrials_pressed(a,:)=smooth(all_mtrials_pressed(a,:),4);
    
%     maxval=max([all_mtrials_missed(a,:) all_mtrials_pressed(a,:)]);
%     all_mtrials_missed(a,:)=all_mtrials_missed(a,:)./maxval;
%     all_mtrials_pressed(a,:)=all_mtrials_pressed(a,:)./maxval;
end

figure;
subplot(1,2,1);
imagesc(all_mtrials_missed(:,t_0:t_5));
colormap(parula);
colorbar
box off;
caxis([0 1]);
subplot(1,2,2);
imagesc(all_mtrials_pressed(:,t_0:t_5));
colormap(parula);
colorbar
box off;
caxis([0 1]);

figure;
imagesc(all_mtrials_missed(:,t_0:t_5)-all_mtrials_pressed(:,t_0:t_5));
colormap(redblue);
colorbar;
maxval=max(max(abs(all_mtrials_missed(:,t_0:t_5)-all_mtrials_pressed(:,t_0:t_5))));
caxis([-maxval maxval])

%ranksum(mean(all_mtrials_missed(:,t_0:t_5),2),mean(all_mtrials_pressed(:,t_0:t_5),2))
shock_miss=mean(all_mtrials_missed(:,t_0:t_5),2);
shock_press=mean(all_mtrials_pressed(:,t_0:t_5),2);
figure; bar([1 2],[mean(shock_miss) mean(shock_press)]);
hold on
errorbar([1 2],[mean(shock_miss) mean(shock_press)],[std(shock_miss) std(shock_press)]./sqrt(length(shock_miss)));

all_mtrials_pressed=[];
all_mtrials_missed=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/results/',filename,'.mat'));
        
    temp=squeeze(mean(dfof_trials_pressed_tone,2));
    all_mtrials_pressed=[all_mtrials_pressed; temp(reward_only,:)];
    temp=squeeze(mean(dfof_trials_missed_tone,2));
    all_mtrials_missed=[all_mtrials_missed; temp(reward_only,:)];

end
t_0=find(time_trial>=0,1,'first');
t_1=find(time_trial>=1,1,'first');
t_2=find(time_trial>=2,1,'first');
t_3=find(time_trial>=3,1,'first');
t_4=find(time_trial>=4,1,'first');
t_5=find(time_trial<5,1,'last');

num_cells=size(all_mtrials_missed,1);

ind=[];
for a=1:num_cells
    all_mtrials_missed(a,:)=smooth(all_mtrials_missed(a,:),4);
    all_mtrials_pressed(a,:)=smooth(all_mtrials_pressed(a,:),4);
    
%     maxval=max([all_mtrials_missed(a,:) all_mtrials_pressed(a,:)]);
%     all_mtrials_missed(a,:)=all_mtrials_missed(a,:)./maxval;
%     all_mtrials_pressed(a,:)=all_mtrials_pressed(a,:)./maxval;
end

figure;
subplot(1,2,1);
imagesc(all_mtrials_missed(:,t_0:t_5));
colormap(parula);
colorbar
box off;
caxis([0 1]);
subplot(1,2,2);
imagesc(all_mtrials_pressed(:,t_0:t_5));
colormap(parula);
colorbar
box off;
caxis([0 1]);

figure;
imagesc(all_mtrials_missed(:,t_0:t_5)-all_mtrials_pressed(:,t_0:t_5));
colormap(redblue);
colorbar;
maxval=max(max(abs(all_mtrials_missed(:,t_0:t_5)-all_mtrials_pressed(:,t_0:t_5))));
caxis([-maxval maxval])

%ranksum(mean(all_mtrials_missed(:,t_0:t_5),2),mean(all_mtrials_pressed(:,t_0:t_5),2))

reward_miss=mean(all_mtrials_missed(:,t_0:t_5),2);
reward_press=mean(all_mtrials_pressed(:,t_0:t_5),2);

figure; bar([1 2],[mean(reward_miss) mean(reward_press)]);
hold on
errorbar([1 2],[mean(reward_miss) mean(reward_press)],[std(reward_miss) std(reward_press)]./sqrt(length(reward_miss)));
