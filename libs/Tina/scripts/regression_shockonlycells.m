function regression_shockonlycells(pathname,filedate,fn)

currpath=strcat(pathname,filedate);

load(strcat(currpath,'/regression/',fn,'_shock.mat'));
shock=behav_cells;
load(strcat(currpath,'/regression/',fn,'_reward.mat'));
reward=behav_cells;
load(strcat(currpath,'/regression/',fn,'_lever.mat'));
lever=behav_cells;

% find shock+reward+lever cells
lever_shock_reward=intersect(intersect(shock,lever),reward);
num_lever_shock_reward=length(lever_shock_reward);
display(num_lever_shock_reward);

% find shock+reward cells
shock_reward=intersect(shock,reward);
shock_reward=shock_reward(~ismember(shock_reward,lever));
num_shock_reward=length(shock_reward);
display(num_shock_reward);

% find lever+shock cells
lever_shock=intersect(shock,lever);
lever_shock=lever_shock(~ismember(lever_shock,reward));
num_lever_shock=length(lever_shock);
display(num_lever_shock);

% find lever+reward cells
lever_reward=intersect(reward,lever);
lever_reward=lever_reward(~ismember(lever_reward,shock));
num_lever_reward=length(lever_reward);
display(num_lever_reward);

% find shock_only cells
overlap=ismember(shock,reward);
del_reward=find(overlap==1);
overlap=ismember(shock,lever);
del_lever=find(overlap==1);
del=unique([del_reward del_lever]);
shock_only=shock;
shock_only(del)=[];

num_shockonly=length(shock_only);
display(num_shockonly);

savename=strcat(currpath,'/regression/',fn,'_shockonly.mat');
save(savename,'shock_only');

% find reward_only cells
overlap=ismember(reward,shock);
del_shock=find(overlap==1);
overlap=ismember(reward,lever);
del_lever=find(overlap==1);
del=unique([del_shock del_lever]);
reward_only=reward;
reward_only(del)=[];

num_rewardonly=length(reward_only);
display(num_rewardonly);

savename=strcat(currpath,'/regression/',fn,'_rewardonly.mat');
save(savename,'reward_only');

% find lever_only cells
overlap=ismember(lever,shock);
del_shock=find(overlap==1);
overlap=ismember(lever,reward);
del_reward=find(overlap==1);
del=unique([del_shock del_reward]);
lever_only=lever;
lever_only(del)=[];

num_leveronly=length(lever_only);
display(num_leveronly);

savename=strcat(currpath,'/regression/',fn,'_leveronly.mat');
save(savename,'lever_only');

load(strcat(currpath,'/timecourses/',fn,'_cellmasks.mat'));
% figure;
% subplot(2,2,1);
% imagesc(movm); 
% hold on
% colormap(gray); 
% for a=1:size(cellmask,1)
%     contour(squeeze(cellmask(a,:,:)),1,'w');
% end
% for a=1:length(shock_only)
%     contour(squeeze(cellmask(shock_only(a),:,:)),1,'r','linewidth',2);
%     hold on
%     text(segcentroid(shock_only(a),1), segcentroid(shock_only(a),2), num2str(shock_only(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% 
% subplot(2,2,2);
% imagesc(movm); 
% hold on
% colormap(gray); 
% for a=1:size(cellmask,1)
%     contour(squeeze(cellmask(a,:,:)),1,'w');
% end
% for a=1:length(reward_only)
%     contour(squeeze(cellmask(reward_only(a),:,:)),1,'c','linewidth',2);
%     hold on
%     text(segcentroid(reward_only(a),1), segcentroid(reward_only(a),2), num2str(reward_only(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% 
% subplot(2,2,3);
% imagesc(movm); 
% hold on
% colormap(gray);
% for a=1:length(reward_only)
%     contour(squeeze(cellmask(reward_only(a),:,:)),1,'c','linewidth',2);
% end
% for a=1:length(shock_only)
%     contour(squeeze(cellmask(shock_only(a),:,:)),1,'r','linewidth',2);
% end
% 
% 
% subplot(2,2,4);
% imagesc(movm); 
% hold on
% colormap(gray); 
% for a=1:size(cellmask,1)
%     contour(squeeze(cellmask(a,:,:)),1,'w');
% end
% for a=1:length(lever_only)
%     contour(squeeze(cellmask(lever_only(a),:,:)),1,'g','linewidth',2);
%     hold on
%     text(segcentroid(lever_only(a),1), segcentroid(lever_only(a),2), num2str(lever_only(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% 
% mean(segcentroid(shock_only,2))
% mean(segcentroid(reward_only,2))


figure;
imagesc(movm); 
hold on
colormap(gray);
for a=1:size(cellmask,1)
    contour(squeeze(cellmask(a,:,:)),1,'w');
end
for a=1:length(reward_only)
    contour(squeeze(cellmask(reward_only(a),:,:)),1,'c','linewidth',2);
%     text(segcentroid(reward_only(a),1), segcentroid(reward_only(a),2), num2str(reward_only(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w');
end
figure;
imagesc(movm); 
hold on
colormap(gray);
for a=1:size(cellmask,1)
    contour(squeeze(cellmask(a,:,:)),1,'w');
end
for a=1:length(shock_only)
    contour(squeeze(cellmask(shock_only(a),:,:)),1,'r','linewidth',2);
end

end