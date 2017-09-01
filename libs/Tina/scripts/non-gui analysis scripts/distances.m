close all
clear all 

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% calculate distances between the shock neurons or reward neurons, and
% ranodmly chosen neurons

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};


% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};


dist_shock_all=[];
mice=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    
    % change: _rewardonly or _shockonly
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    
    % change: _rewardonly or _shockonly
    shock_centroids=segcentroid(reward_only,:);
    
    dist_shock=dist(shock_centroids');
    % distances in m/l direction
%     dist_shock=dist(shock_centroids(:,2)');
    dist_shock=triu(dist_shock,1);
    dist_shock=dist_shock(dist_shock~=0);
    dist_shock_all=[dist_shock_all; dist_shock];
    mice=[mice mean(dist_shock)];
   
    % change: _rewardonly or _shockonly
    choose_shock(z)=length(reward_only);
end

m_dist_shock=mean(dist_shock_all);

nshuffs=1000;
m_dist_nonmod=zeros(1,nshuffs);
for p=1:nshuffs
    dist_nonmod_all=[];
    for z=1:length(dates)
        date=dates{z};
        filename=filenames{z};
        
        load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
        nonmod_centroids=segcentroid;
        choose_nonmod=randperm(length(nonmod_centroids),choose_shock(z));
        nonmod_centroids=nonmod_centroids(choose_nonmod,:);
        
        dist_nonmod=dist(nonmod_centroids');
        % distances in m/l direction
%         dist_nonmod=dist(nonmod_centroids(:,2)');
        dist_nonmod=triu(dist_nonmod,2);
        dist_nonmod=dist_nonmod(dist_nonmod~=0);
        dist_nonmod_all=[dist_nonmod_all; dist_nonmod];
    end
    m_dist_nonmod(p)=mean(dist_nonmod_all);
    p
end

pval=length(find(m_dist_nonmod<m_dist_shock))/nshuffs

figure;
plot([1,2],[m_dist_shock,mean(m_dist_nonmod)],'.-')
hold on
errorbar([1,2],[m_dist_shock,mean(m_dist_nonmod)],[std(dist_shock_all)/sqrt(length(dist_shock_all)),std(m_dist_nonmod)/sqrt(length(m_dist_nonmod))]);

% %%
% close all
% clear all
% pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
% 
% % calculate distances between the shock neurons, reward neurons, and
% % non-mod neurons
% 
% % dates={'20160429','20160429','20160429','20160429','20160528'};
% % filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};
% 
% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};
% 
% 
% dist_shock_all=[];
% dist_shockneg_all=[];
% dist_reward_all=[];
% dist_rewardneg_all=[];
% for z=1:length(dates)
%     date=dates{z};
%     filename=filenames{z};
% 
%     load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
%     load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
% 
%     shock_centroids=segcentroid(shock_only,:);
%     dist_shock=dist(shock_centroids');
%     dist_shock=triu(dist_shock,1);
%     dist_shock=dist_shock(dist_shock~=0);
%     dist_shock_all=[dist_shock_all; dist_shock];
%     
%     load(strcat(pathname,date,'/regression/',filename,'_shockonly_neg.mat'));
% 
%     shockneg_centroids=segcentroid(shock_only,:);
%     dist_shockneg=dist(shockneg_centroids');
%     dist_shockneg=triu(dist_shockneg,1);
%     dist_shockneg=dist_shockneg(dist_shockneg~=0);
%     dist_shockneg_all=[dist_shockneg_all; dist_shockneg];
%     
%     load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
% 
%     reward_centroids=segcentroid(reward_only,:);
%     dist_reward=dist(reward_centroids');
%     dist_reward=triu(dist_reward,1);
%     dist_reward=dist_reward(dist_reward~=0);
%     dist_reward_all=[dist_reward_all; dist_reward];
%     
%     load(strcat(pathname,date,'/regression/',filename,'_rewardonly_neg.mat'));
%     
%     rewardneg_centroids=segcentroid(reward_only_neg,:);
%     dist_rewardneg=dist(rewardneg_centroids');
%     dist_rewardneg=triu(dist_rewardneg,1);
%     dist_rewardneg=dist_rewardneg(dist_rewardneg~=0);
%     dist_rewardneg_all=[dist_rewardneg_all; dist_rewardneg];
% 
% end
% m_dist_shock=mean(dist_shock_all);
% m_dist_shockneg=mean(dist_shockneg_all);
% m_dist_reward=mean(dist_reward_all);
% m_dist_rewardneg=mean(dist_rewardneg_all);
% 
% sem_dist_shock=std(dist_shock_all)/sqrt(length(dist_shock_all));
% sem_dist_shockneg=std(dist_shockneg_all)/sqrt(length(dist_shockneg_all));
% sem_dist_reward=std(dist_reward_all)/sqrt(length(dist_reward_all));
% sem_dist_rewardneg=std(dist_rewardneg_all)/sqrt(length(dist_rewardneg_all));
% 
% figure;
% bar([1:4],[m_dist_shock,m_dist_reward,m_dist_shockneg,m_dist_rewardneg]);
% hold on
% errorbar([1:4],[m_dist_shock,m_dist_reward,m_dist_shockneg,m_dist_rewardneg],[sem_dist_shock,sem_dist_reward,sem_dist_shockneg,sem_dist_rewardneg]);
% 
% [~,~,stats]=anova1([dist_shock_all' dist_reward_all' dist_shockneg_all' dist_rewardneg_all'],[zeros(1,length(dist_shock_all)) ones(1,length(dist_reward_all)) ones(1,length(dist_shockneg_all))*2 ones(1,length(dist_rewardneg_all))*3]);
% [c,~,~,gnames] = multcompare(stats);