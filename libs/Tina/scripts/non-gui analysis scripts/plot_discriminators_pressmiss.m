close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

all_heatmap_pressed_pressed=[];
all_heatmap_pressed_missed=[];
all_heatmap_pressed_shock=[];
all_heatmap_pressed_reward=[];

all_heatmap_missed_pressed=[];
all_heatmap_missed_missed=[];
all_heatmap_missed_shock=[];
all_heatmap_missed_reward=[];

total_cells_pressed=[];
total_cells_missed=[];

for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/results/',filename,'.mat'));
    load(strcat(pathname,date,'/results/',filename,'_pressmiss.mat'));

    
    t_0=find(time_trial>=0,1,'first');
    t_5=find(time_trial<5,1,'last');
    
    num_cells=length(discriminators);
    
    diff_discriminators=mean(mtrials_pressed_tone(discriminators,t_0:t_5)-mtrials_missed_tone(discriminators,t_0:t_5),2);
    
    missed_cells=find(diff_discriminators<0);
    pressed_cells=find(diff_discriminators>0);

    mtrials_missed_tone=mtrials_missed_tone(discriminators,:);
    mtrials_pressed_tone=mtrials_pressed_tone(discriminators,:);
    mtrials_shock_press=mtrials_shock_press(discriminators,:);
    mtrials_reward_press=mtrials_reward_press(discriminators,:);
    
    % normalize each cell by max value
    for a=1:size(mtrials_missed_tone,1)
        maxval_missed=max(mtrials_missed_tone(a,t_0:t_5));
        maxval_pressed=max(mtrials_pressed_tone(a,t_0:t_5));
        maxval=max([maxval_missed maxval_pressed]);
        mtrials_missed_tone(a,:)=mtrials_missed_tone(a,:)./maxval;
        mtrials_pressed_tone(a,:)=mtrials_pressed_tone(a,:)./maxval;
        maxval_reward=max(mtrials_reward_press(a,t_0:t_5));
        maxval_shock=max(mtrials_shock_press(a,t_0:t_5));
        maxval=max([maxval_reward maxval_shock]);
        mtrials_shock_press(a,:)=mtrials_shock_press(a,:)./maxval;
        mtrials_reward_press(a,:)=mtrials_reward_press(a,:)./maxval;
    end
        
    heatmap_pressed_pressed=mtrials_pressed_tone(pressed_cells,t_0:t_5);
    heatmap_pressed_missed=mtrials_missed_tone(pressed_cells,t_0:t_5);
    heatmap_pressed_reward=mtrials_reward_press(pressed_cells,t_0:t_5);
    heatmap_pressed_shock=mtrials_shock_press(pressed_cells,t_0:t_5);
    
    heatmap_missed_pressed=mtrials_pressed_tone(missed_cells,t_0:t_5);
    heatmap_missed_missed=mtrials_missed_tone(missed_cells,t_0:t_5);
    heatmap_missed_reward=mtrials_reward_press(missed_cells,t_0:t_5);
    heatmap_missed_shock=mtrials_shock_press(missed_cells,t_0:t_5);
    
    all_heatmap_pressed_pressed=[all_heatmap_pressed_pressed; heatmap_pressed_pressed];
    all_heatmap_pressed_missed=[all_heatmap_pressed_missed; heatmap_pressed_missed];
    all_heatmap_pressed_reward=[all_heatmap_pressed_reward; heatmap_pressed_reward];
    all_heatmap_pressed_shock=[all_heatmap_pressed_shock; heatmap_pressed_shock];
    
    all_heatmap_missed_pressed=[all_heatmap_missed_pressed; heatmap_missed_pressed];
    all_heatmap_missed_missed=[all_heatmap_missed_missed; heatmap_missed_missed];
    all_heatmap_missed_reward=[all_heatmap_missed_reward; heatmap_missed_reward];
    all_heatmap_missed_shock=[all_heatmap_missed_shock; heatmap_missed_shock];
    
    total_cells_pressed(z)=length(pressed_cells);
    total_cells_missed(z)=length(missed_cells);
end
time_axis=time_trial(t_0:t_5);
t_1=find(time_axis>=1,1,'first');    

border_pressed(1)=total_cells_pressed(1);
border_missed(1)=total_cells_missed(1);
for a=2:length(dates)
    border_pressed(a)=total_cells_pressed(a)+border_pressed(a-1);
    border_missed(a)=total_cells_missed(a)+border_missed(a-1);
end
border_pressed=border_pressed+0.5;
border_missed=border_missed+0.5;

figure;
subplot(1,4,1);
imagesc((all_heatmap_missed_missed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
subplot(1,4,2);
imagesc((all_heatmap_missed_pressed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
subplot(1,4,3);
imagesc((all_heatmap_missed_shock));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
plot([t_1 t_1],[0.5 size(all_heatmap_missed_shock,1)+0.5],'-w','linewidth',2);
subplot(1,4,4);
imagesc((all_heatmap_missed_reward));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
plot([t_1 t_1],[0.5 size(all_heatmap_missed_shock,1)+0.5],'-w','linewidth',2);


figure;
subplot(1,4,1);
imagesc((all_heatmap_pressed_missed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_pressed' border_pressed'],'-w','linewidth',2);
subplot(1,4,2);
imagesc((all_heatmap_pressed_pressed));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_pressed' border_pressed'],'-w','linewidth',2);
subplot(1,4,3);
imagesc((all_heatmap_pressed_shock));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
plot([t_1 t_1],[0.5 size(all_heatmap_pressed_shock,1)+0.5],'-w','linewidth',2);
subplot(1,4,4);
imagesc((all_heatmap_pressed_reward));
colorbar
caxis([0 1])
hold on
plot([0.5 38.5],[border_missed' border_missed'],'-w','linewidth',2);
plot([t_1 t_1],[0.5 size(all_heatmap_pressed_shock,1)+0.5],'-w','linewidth',2);

