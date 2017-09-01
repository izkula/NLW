%close all
clear all 

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% calculate position of shock neurons or reward neurons, and
% ranodmly chosen neurons

% dates={'20160429','20160429','20160429','20160429','20160528'};
% filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

dates={'20160529','20160529','20160730','20160730','20160731'};
filenames={'GRIN_ pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};


centroid_shock_all=[];
mice=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    % change: _rewardonly or _shockonly
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));

    % change: reward_only or shock_only
    centroid_shock=segcentroid(reward_only,1);
    centroid_shock_all=[centroid_shock_all; centroid_shock];
    mice=[mice mean(centroid_shock)];
    % change: reward_only or shock_only
    choose_shock(z)=length(reward_only);
end

m_centroid_shock=mean(centroid_shock_all);

numshuffs=1000;
m_centroid_nonmod=zeros(1,numshuffs);
for p=1:numshuffs
    centroid_nonmod_all=[];
    for z=1:length(dates)
        date=dates{z};
        filename=filenames{z};
        
        load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
        centroid_nonmod=segcentroid(:,1);
        choose_nonmod=randperm(length(centroid_nonmod),choose_shock(z));
        centroid_nonmod=centroid_nonmod(choose_nonmod);
        centroid_nonmod_all=[centroid_nonmod_all; centroid_nonmod];
    end
    m_centroid_nonmod(p)=mean(centroid_nonmod_all);
    p
end

pval=length(find(m_centroid_nonmod<m_centroid_shock))/numshuffs
pval=length(find(m_centroid_nonmod>m_centroid_shock))/numshuffs

figure;
plot([1],[m_centroid_shock],'.-k','markersize',20)
hold on
errorbar([1],[m_centroid_shock],[std(centroid_shock_all)/sqrt(length(centroid_shock_all))],'k');
plot([1],[m_centroid_nonmod],'.r')

figure;
histogram(m_centroid_nonmod,'binwidth',2.5);
hold on
plot([m_centroid_shock m_centroid_shock],[0 200],'-r');

% figure;
% plot([1,2],[m_centroid_shock,mean(m_centroid_nonmod)],'.-')
% hold on
% errorbar([1,2],[m_centroid_shock,mean(m_centroid_nonmod)],[std(centroid_shock_all)/sqrt(length(centroid_shock_all)),std(m_centroid_nonmod)/sqrt(length(m_centroid_nonmod))]);

%%
close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% calculate position of shock neurons or reward neurons, and
% ranodmly chosen neurons

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};


centroid_shock_all=[];
% centroid_shockneg_all=[];
centroid_reward_all=[];
% centroid_rewardneg_all=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));

    centroid_shock=segcentroid(shock_only,1);
    centroid_shock_all=[centroid_shock_all; centroid_shock];
    
%     load(strcat(pathname,date,'/regression/',filename,'_shockonly_neg.mat'));
% 
%     centroid_shockneg=segcentroid(shock_only_neg,1);
%     centroid_shockneg_all=[centroid_shockneg_all; centroid_shockneg];
    
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));

    centroid_reward=segcentroid(reward_only,1);
    centroid_reward_all=[centroid_reward_all; centroid_reward];
    
%     load(strcat(pathname,date,'/regression/',filename,'_rewardonly_neg.mat'));
% 
%     centroid_rewardneg=segcentroid(reward_only_neg,1);
%     centroid_rewardneg_all=[centroid_rewardneg_all; centroid_rewardneg];

    mice_shock(z)=mean(centroid_shock);
    mice_reward(z)=mean(centroid_reward);
    
end
m_centroid_shock=mean(centroid_shock_all);
% m_centroid_shockneg=mean(centroid_shockneg_all);
m_centroid_reward=mean(centroid_reward_all);
% m_centroid_rewardneg=mean(centroid_rewardneg_all);

sem_centroid_shock=std(centroid_shock_all)/sqrt(length(centroid_shock_all));
% sem_centroid_shockneg=std(centroid_shockneg_all)/sqrt(length(centroid_shockneg_all));
sem_centroid_reward=std(centroid_reward_all)/sqrt(length(centroid_reward_all));
% sem_centroid_rewardneg=std(centroid_rewardneg_all)/sqrt(length(centroid_rewardneg_all));


figure;
bar([1:2],[m_centroid_shock,m_centroid_reward]);
hold on
errorbar([1:2],[m_centroid_shock,m_centroid_reward],[sem_centroid_shock,sem_centroid_reward]);

[~,~,stats]=anova1([centroid_shock_all' centroid_reward_all'],[zeros(1,length(centroid_shock_all)) ones(1,length(centroid_reward_all))]);
[c,~,~,gnames] = multcompare(stats);

figure;
plot([1 2],[mice_shock' mice_reward'],'.-','markersize',30);


% figure;
% bar([1:4],[m_centroid_shock,m_centroid_reward,m_centroid_shockneg,m_centroid_rewardneg]);
% hold on
% errorbar([1:4],[m_centroid_shock,m_centroid_reward,m_centroid_shockneg,m_centroid_rewardneg],[sem_centroid_shock,sem_centroid_reward,sem_centroid_shockneg,sem_centroid_rewardneg]);
% 
% [~,~,stats]=anova1([centroid_shock_all' centroid_reward_all' centroid_shockneg_all' centroid_rewardneg_all'],[zeros(1,length(centroid_shock_all)) ones(1,length(centroid_reward_all)) ones(1,length(centroid_shockneg_all))*2 ones(1,length(centroid_rewardneg_all))*3]);
% [c,~,~,gnames] = multcompare(stats);