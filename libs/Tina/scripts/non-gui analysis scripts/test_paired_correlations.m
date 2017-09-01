close all
clear all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};


allshockcorr=[];
allshockrewardcorr=[];
allshocknonmodcorr=[];
allrewardcorr=[];
allrewardnonmodcorr=[];
allnonmodcorr=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/timecourses/',filename,'_dfof.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    load(strcat(pathname,date,'/behavior/',filename,'.mat'));
    
    pre_ind=find(time>=2,1,'first');
    post_ind=find(time>=10,1,'first');
    
    del_trials=[];
    for a=1:length(tone_start)
        tone_start_2p(a)=find(time<=time_behav(tone_start(a)),1,'last')-pre_ind;
        tone_end_2p(a)=tone_start_2p(a)+post_ind;
        del_trials=[del_trials tone_start_2p(a):tone_end_2p(a)];
    end
    
    dfof_sub(:,del_trials)=[];
    allcorr=corr(dfof_sub');
    
    % get non-modulated neurons
    load(strcat(pathname,date,'/regression/',filename,'_lever.mat'));
    all_lever=behav_cells;
    load(strcat(pathname,date,'/regression/',filename,'_reward.mat'));
    all_reward=behav_cells;
    load(strcat(pathname,date,'/regression/',filename,'_shock.mat'));
    all_shock=behav_cells;
        
    task_cells=unique([all_lever all_reward all_shock]);
    nonmod=[1:size(dfof_sub),1];
    nonmod(task_cells)=[];
    
    % find shock-shock correlations
    shockcorr=corr(dfof_sub(shock_only,:)');
    shockcorr=nonzeros(triu(shockcorr,1)');
    %allshockcorr=[allshockcorr; shockcorr];
    allshockcorr=[allshockcorr; mean(shockcorr)];
    
    % find shock-to-nonshock
    shocknonmodcorr=[];
    for a=1:length(shock_only)
        curr_cell=shock_only(a);      
        shocknonmodcorr=[shocknonmodcorr; allcorr(curr_cell,nonmod)'];
    end
    %allshocknonmodcorr=[allshocknonmodcorr; shocknonmodcorr];
    allshocknonmodcorr=[allshocknonmodcorr; mean(shocknonmodcorr)];

    if ~isempty(reward_only)
        % find correlations between reward-reward neurons
        rewardcorr=corr(dfof_sub(reward_only,:)');
        rewardcorr=nonzeros(triu(rewardcorr,1)');
        %allrewardcorr=[allrewardcorr; rewardcorr];
        allrewardcorr=[allrewardcorr; mean(rewardcorr)];
        
        % find correlations between reward-nonmod neurons
        rewardnonmodcorr=[];
        for a=1:length(reward_only)
            curr_cell=reward_only(a);
            rewardnonmodcorr=[rewardnonmodcorr; allcorr(curr_cell,nonmod)'];
        end
        %allrewardnonmodcorr=[allrewardnonmodcorr; rewardnonmodcorr];
        allrewardnonmodcorr=[allrewardnonmodcorr; mean(rewardnonmodcorr)];
        
        % find correlations between reward-shock neurons
        shockrewardcorr=[];
        for a=1:length(shock_only)
            curr_cell=shock_only(a);
            shockrewardcorr=[shockrewardcorr; allcorr(curr_cell,reward_only)'];
        end
        %allshockrewardcorr=[allshockrewardcorr; shockrewardcorr];
        allshockrewardcorr=[allshockrewardcorr; mean(shockrewardcorr)];
    end
    
    % find correlations between nonmod-nonmod neurons
    nonmodcorr=corr(dfof_sub(nonmod,:)');
    nonmodcorr=nonzeros(triu(nonmodcorr,1)');
    %allnonmodcorr=[allnonmodcorr; nonmodcorr];
    allnonmodcorr=[allnonmodcorr; mean(nonmodcorr)];
end


% [~,~,stats]=anova1([allshockcorr' allrewardcorr' allnonmodcorr'],[zeros(1,length(allshockcorr)) ones(1,length(allrewardcorr)) ones(1,length(allnonmodcorr))*2]);
% [c,~,~,gnames] = multcompare(stats);
% 
mshock=mean(allshockcorr);
mreward=mean(allrewardcorr);
mnonmod=mean(allnonmodcorr);
mshocknonmod=mean(allshocknonmodcorr);
mrewardnonmod=mean(allrewardnonmodcorr);
mshockreward=mean(allshockrewardcorr);

semshock=std(allshockcorr)/sqrt(length(allshockcorr));
semreward=std(allrewardcorr)/sqrt(length(allrewardcorr));
semnonmod=std(allnonmodcorr)/sqrt(length(allnonmodcorr));
semshocknonmod=std(allshocknonmodcorr)/sqrt(length(allshocknonmodcorr));
semrewardnonmod=std(allrewardnonmodcorr)/sqrt(length(allrewardnonmodcorr));
semshockreward=std(allshockrewardcorr)/sqrt(length(allshockrewardcorr));

figure;
bar([1:3],[mshock,mshocknonmod,mshockreward]);
hold on
errorbar([1:3],[mshock,mshocknonmod,mshockreward],[semshock,semshocknonmod,semshockreward]);

figure;
bar([1:3],[mreward,mrewardnonmod,mshockreward]);
hold on
errorbar([1:3],[mreward,mrewardnonmod,mshockreward],[semreward,semrewardnonmod,semshockreward]);

% [~,~,stats]=anova1([allshockcorr' allshocknonmodcorr' allnonmodcorr'],[zeros(1,length(allshockcorr)) ones(1,length(allrewardcorr)) ones(1,length(allnonmodcorr))*2]);
% [c,~,~,gnames] = multcompare(stats);

% %%
% 
% close all
% clear all
% 
% pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
% 
% dates={'20160429','20160429','20160429','20160429','20160528'};
% filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};
% 
% 
% allshockcorr=[];
% allrewardcorr=[];
% allnonmodcorr=[];
% 
% for z=1:length(dates)
% %for z=1:4
%     date=dates{z};
%     filename=filenames{z};
%     
%     load(strcat(pathname,date,'/timecourses/',filename,'_dfof.mat'));
%     load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
%     load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
%     load(strcat(pathname,date,'/behavior/',filename,'.mat'));
%     
%     pre_ind=find(time>=2,1,'first');
%     post_ind=find(time>=10,1,'first');
%     
%     del_trials=[];
%     for a=1:length(tone_start)
%         tone_start_2p(a)=find(time<=time_behav(tone_start(a)),1,'last')-pre_ind;
%         tone_end_2p(a)=tone_start_2p(a)+post_ind;
%         del_trials=[del_trials tone_start_2p(a):tone_end_2p(a)];
%     end
%     
%     dfof_sub(:,del_trials)=[];
%     allcorr=corr(dfof_sub');
%     
%     % load all cells modulated (even shock+lever etc)
%     load(strcat(pathname,date,'/regression/',filename,'_lever.mat'));
%     all_lever=behav_cells;
%     load(strcat(pathname,date,'/regression/',filename,'_reward.mat'));
%     all_reward=behav_cells;
%     load(strcat(pathname,date,'/regression/',filename,'_shock.mat'));
%     all_shock=behav_cells;
%         
%     task_cells=unique([all_lever all_reward all_shock]);
%     
%     nonmod=[1:size(dfof_sub),1];
%     nonmod(task_cells)=[];
%     
% %     % find shock-to-nonshock
% %     nonshockcorr=[];
% %     for a=1:length(shock_only)
% %         curr_cell=shock_only(a);
% %         
% %         other_shock=shock_only;
% %         other_shock(a)=[];
% %         
% %         non_shock=[1:size(dfof_sub,1)];
% %         non_shock(shock_only)=[];
% %         
% %         nonshockcorr=[nonshockcorr allcorr(curr_cell,non_shock)];
% %     end
%     
%     % find correlations between shock-shock neurons
%     shockcorr=corr(dfof_sub(shock_only,:)');
%     %shockcorr=triu(shockcorr,1);
%     shockcorr=nonzeros(triu(shockcorr,1)');
%     allshockcorr=[allshockcorr; shockcorr];
% 
%     
%     if ~isempty(reward_only)
%         % find correlations between reward-reward neurons
%         rewardcorr=corr(dfof_sub(reward_only,:)');
%         %shockcorr=triu(shockcorr,1);
%         rewardcorr=nonzeros(triu(rewardcorr,1)');
%         
%         allrewardcorr=[allrewardcorr; rewardcorr];
%     end
%     
%     % find correlations between nonmod-nonmod neurons
%     nonmodcorr=corr(dfof_sub(nonmod,:)');
%     %shockcorr=triu(shockcorr,1);
%     nonmodcorr=nonzeros(triu(nonmodcorr,1)');
%     
%     allnonmodcorr=[allnonmodcorr; nonmodcorr];
%     %allnonshockcorr=[allnonshockcorr nonshockcorr];
% end
% 
% [~,~,stats]=anova1([allshockcorr' allrewardcorr' allnonmodcorr'],[zeros(1,length(allshockcorr)) ones(1,length(allrewardcorr)) ones(1,length(allnonmodcorr))*2]);
% [c,~,~,gnames] = multcompare(stats);
% 
% mshock=mean(allshockcorr);
% mreward=mean(allrewardcorr);
% mnonmod=mean(allnonmodcorr);
% 
% semshock=std(allshockcorr)/sqrt(length(allshockcorr));
% semreward=std(allrewardcorr)/sqrt(length(allrewardcorr));
% semnonmod=std(allnonmodcorr)/sqrt(length(allnonmodcorr));
% 
% figure;
% bar([1:3],[mshock,mreward,mnonmod]);
% hold on
% errorbar([1:3],[mshock,mreward,mnonmod],[semshock,semreward,semnonmod]);
% 
