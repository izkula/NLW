clear all
close all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';


% dates={'20160429','20160429','20160429','20160429','20160528'};
% filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

dates={'20160529','20160529','20160730','20160730','20160731'};
filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

all_mtrials_reward=[];
all_mtrials_shock=[];
all_mtrials_lever=[];

for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_leveronly.mat'));
    load(strcat(pathname,date,'/results/',filename,'.mat'));
    
    mtrials_lever_press=[mtrials_reward_press;mtrials_shock_press];
    
    all_mtrials_shock=[all_mtrials_shock; mtrials_shock_press(shock_only,:)];
    all_mtrials_reward=[all_mtrials_reward; mtrials_reward_press(reward_only,:)];
    all_mtrials_lever=[all_mtrials_lever; mtrials_lever_press(lever_only,:)];

end

%%
t_0=find(time_trial>=0,1,'first');
t_1=find(time_trial>=1,1,'first');
t_2=find(time_trial>=2,1,'first');
t_3=find(time_trial>=3,1,'first');
t_4=find(time_trial>=4,1,'first');

num_shock=size(all_mtrials_shock,1);
num_reward=size(all_mtrials_reward,1);
num_lever=size(all_mtrials_lever,1);

%%

ind=[];
for a=1:num_shock
    all_mtrials_shock(a,:)=smooth(all_mtrials_shock(a,:),4);
%     all_mtrials_shock(a,:)=all_mtrials_shock(a,:)./max(all_mtrials_shock(a,:));
    [~,ind(a)]=max(all_mtrials_shock(a,t_1:t_3));
end

order=[ind' [1:num_shock]'];
order=sortrows(order,1);

all_mtrials_shock_sort=[];
for a=1:length(order)
    all_mtrials_shock_sort(a,:)=all_mtrials_shock(order(a,2),:);
end


figure;
imagesc(flipud(all_mtrials_shock_sort));
colormap(parula);
colorbar
box off;
caxis([0 1.5]);
hold on
plot([t_1 t_1],[0 num_shock+1],'-w','linewidth',2)
plot([t_0 t_0],[0 num_shock+1],'-w','linewidth',2)

%%
ind=[];
for a=1:num_reward
    all_mtrials_reward(a,:)=smooth(all_mtrials_reward(a,:),4);
%     all_mtrials_reward(a,:)=all_mtrials_reward(a,:)./max(all_mtrials_reward(a,:));
    [~,ind(a)]=max(all_mtrials_reward(a,t_1:t_3));
end

% order=[ind' [1:num_reward]'];
% order=sortrows(order,1);

all_mtrials_reward_sort=[];
for a=1:length(order)
    all_mtrials_reward_sort(a,:)=all_mtrials_reward(order(a,2),:);
end


figure;
imagesc(flipud(all_mtrials_reward_sort));
colormap(parula);
colorbar
box off;
caxis([0 1.5]);
hold on
plot([t_1 t_1],[0 num_reward+1],'-w','linewidth',2)
plot([t_0 t_0],[0 num_reward+1],'-w','linewidth',2)

% %%
% ind=[];
% for a=1:num_lever
%     all_mtrials_lever(a,:)=smooth(all_mtrials_lever(a,:),4);
%     all_mtrials_lever(a,:)=all_mtrials_lever(a,:)./max(all_mtrials_lever(a,:));
%     [~,ind(a)]=max(all_mtrials_lever(a,:));
% end
% 
% order=[ind' [1:num_lever]'];
% order=sortrows(order,1);
% 
% all_mtrials_lever_sort=[];
% for a=1:length(order)
%     all_mtrials_lever_sort(a,:)=all_mtrials_lever(order(a,2),:);
% end
% 
% 
% figure;
% imagesc(flipud(all_mtrials_lever_sort));
% colormap(jet);
% colorbar
% box off;
% caxis([0 1]);
% hold on
% plot([t_1 t_1],[0 num_lever+1],'-w','linewidth',2)
% plot([t_0 t_0],[0 num_lever+1],'-w','linewidth',2)
% 
