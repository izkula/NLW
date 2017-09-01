function regression_shockonlycells_neg(pathname,filedate,fn)

display('includ analysis of NEGATIVE correlations');

currpath=strcat(pathname,filedate);

% load data
load(strcat(currpath,'/regression/',fn,'_shock.mat'));
shock=behav_cells;
load(strcat(currpath,'/regression/',fn,'_reward.mat'));
reward=behav_cells;
load(strcat(currpath,'/regression/',fn,'_lever.mat'));
lever=behav_cells;
load(strcat(currpath,'/regression/',fn,'_shock_neg.mat'));
shockneg=behav_cells;
load(strcat(currpath,'/regression/',fn,'_reward_neg.mat'));
rewardneg=behav_cells;
load(strcat(currpath,'/regression/',fn,'_lever_neg.mat'));
leverneg=behav_cells;


% find shock_only cells
overlap=ismember(shock,reward);
del_reward=find(overlap==1);
overlap=ismember(shock,lever);
del_lever=find(overlap==1);
overlap=ismember(shock,rewardneg);
del_rewardneg=find(overlap==1);
overlap=ismember(shock,leverneg);
del_leverneg=find(overlap==1);
del=unique([del_reward del_lever del_rewardneg del_leverneg]);
shock_only=shock;
shock_only(del)=[];

num_shockonly=length(shock_only);
display(num_shockonly);

% savename=strcat(currpath,'/regression/',fn,'_shockonly.mat');
% save(savename,'shock_only');


% find shock_only neg cells
overlap=ismember(shockneg,reward);
del_reward=find(overlap==1);
overlap=ismember(shockneg,lever);
del_lever=find(overlap==1);
overlap=ismember(shockneg,rewardneg);
del_rewardneg=find(overlap==1);
overlap=ismember(shockneg,leverneg);
del_leverneg=find(overlap==1);
del=unique([del_reward del_lever del_rewardneg del_leverneg]);
shock_only_neg=shockneg;
shock_only_neg(del)=[];

num_shockonly_neg=length(shock_only_neg);
display(num_shockonly_neg);

% savename=strcat(currpath,'/regression/',fn,'_shockonly_neg.mat');
% save(savename,'shock_only_neg');


% find reward_only cells
overlap=ismember(reward,shock);
del_shock=find(overlap==1);
overlap=ismember(reward,lever);
del_lever=find(overlap==1);
overlap=ismember(reward,shockneg);
del_shockneg=find(overlap==1);
overlap=ismember(reward,leverneg);
del_leverneg=find(overlap==1);
del=unique([del_shock del_lever del_shockneg del_leverneg]);
reward_only=reward;
reward_only(del)=[];

num_rewardonly=length(reward_only);
display(num_rewardonly);

% savename=strcat(currpath,'/regression/',fn,'_rewardonly.mat');
% save(savename,'reward_only');

% find reward_only neg cells
overlap=ismember(rewardneg,shock);
del_shock=find(overlap==1);
overlap=ismember(rewardneg,lever);
del_lever=find(overlap==1);
overlap=ismember(rewardneg,shockneg);
del_shockneg=find(overlap==1);
overlap=ismember(rewardneg,leverneg);
del_leverneg=find(overlap==1);
del=unique([del_shock del_lever del_shockneg del_leverneg]);
reward_only_neg=rewardneg;
reward_only_neg(del)=[];

num_rewardonly_neg=length(reward_only_neg);
display(num_rewardonly_neg);

% savename=strcat(currpath,'/regression/',fn,'_rewardonly_neg.mat');
% save(savename,'reward_only_neg');


% find lever_only cells
overlap=ismember(lever,shock);
del_shock=find(overlap==1);
overlap=ismember(lever,reward);
del_reward=find(overlap==1);
overlap=ismember(lever,shockneg);
del_shockneg=find(overlap==1);
overlap=ismember(lever,rewardneg);
del_rewardneg=find(overlap==1);
del=unique([del_shock del_reward del_shockneg del_rewardneg]);
lever_only=lever;
lever_only(del)=[];

num_leveronly=length(lever_only);
display(num_leveronly);

% savename=strcat(currpath,'/regression/',fn,'_leveronly.mat');
% save(savename,'lever_only');

% find lever_only neg cells
overlap=ismember(leverneg,shock);
del_shock=find(overlap==1);
overlap=ismember(leverneg,reward);
del_reward=find(overlap==1);
overlap=ismember(leverneg,shockneg);
del_shockneg=find(overlap==1);
overlap=ismember(leverneg,rewardneg);
del_rewardneg=find(overlap==1);
del=unique([del_shock del_reward del_shockneg del_rewardneg]);
lever_only_neg=leverneg;
lever_only_neg(del)=[];

num_leveronly_neg=length(lever_only_neg);
display(num_leveronly_neg);

% savename=strcat(currpath,'/regression/',fn,'_leveronly_neg.mat');
% save(savename,'lever_only_neg');


% load(strcat(currpath,'/timecourses/',fn,'_cellmasks.mat'));
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
% for a=1:length(shock_only_neg)
%     contour(squeeze(cellmask(shock_only_neg(a),:,:)),1,'m','linewidth',2);
%     hold on
%     text(segcentroid(shock_only_neg(a),1), segcentroid(shock_only_neg(a),2), num2str(shock_only_neg(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
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
%     contour(squeeze(cellmask(reward_only(a),:,:)),1,'b','linewidth',2);
%     hold on
%     text(segcentroid(reward_only(a),1), segcentroid(reward_only(a),2), num2str(reward_only(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% for a=1:length(reward_only_neg)
%     contour(squeeze(cellmask(reward_only_neg(a),:,:)),1,'c','linewidth',2);
%     hold on
%     text(segcentroid(reward_only_neg(a),1), segcentroid(reward_only_neg(a),2), num2str(reward_only_neg(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% 
% subplot(2,2,3);
% imagesc(movm); 
% hold on
% colormap(gray); 
% for a=1:length(reward_only)
%     contour(squeeze(cellmask(reward_only(a),:,:)),1,'b','linewidth',2);
% end
% for a=1:length(shock_only)
%     contour(squeeze(cellmask(shock_only(a),:,:)),1,'r','linewidth',2);
% end
% for a=1:length(reward_only_neg)
%     contour(squeeze(cellmask(reward_only_neg(a),:,:)),1,'c','linewidth',2);
% end
% for a=1:length(shock_only_neg)
%     contour(squeeze(cellmask(shock_only_neg(a),:,:)),1,'m','linewidth',2);
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
% for a=1:length(lever_only_neg)
%     contour(squeeze(cellmask(lever_only_neg(a),:,:)),1,'y','linewidth',2);
%     hold on
%     text(segcentroid(lever_only_neg(a),1), segcentroid(lever_only_neg(a),2), num2str(lever_only_neg(a)), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end


end