%close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% dates={'20160429'};
% filenames={'GRIN_pfcnac4_tonerewardshock_day5'};

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};
 
% dates={'20160730'};
% filenames={'GRIN_pfcvta1_tonerewardshock_day4'};

euc_dist=zeros(length(dates),54);
shock_si=[];
reward_si=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    load(strcat(pathname,date,'/results/',filename,'_pca.mat'));
    
    euc_dist(z,:)=sqrt( (mpca_shock(1,:)-mpca_reward(1,:)).^2 + (mpca_shock(2,:)-mpca_reward(2,:)).^2 + (mpca_shock(3,:)-mpca_reward(3,:)).^2 );

    curr_reward=[];
    ntrials=size(pca_trials_reward,2);
    for a=1:ntrials
        curr_trial=squeeze(pca_trials_reward(:,a,:));
        pca_trials_reward_omit=pca_trials_reward;
        pca_trials_reward_omit(:,a,:)=[];
        mpca_reward_omit=squeeze(mean(pca_trials_reward_omit,2));
        d_opp=sqrt( (mpca_shock(1,:)-curr_trial(1,:)).^2 + (mpca_shock(2,:)-curr_trial(2,:)).^2 + (mpca_shock(3,:)-curr_trial(3,:)).^2 );
        d_same=sqrt( (mpca_reward_omit(1,:)-curr_trial(1,:)).^2 + (mpca_reward_omit(2,:)-curr_trial(2,:)).^2 + (mpca_reward_omit(3,:)-curr_trial(3,:)).^2 );
        curr_reward=[curr_reward; (d_opp-d_same)./(d_opp+d_same)];
    end
    reward_si(z,:)=mean(curr_reward);
    curr_shock=[];
    ntrials=size(pca_trials_shock,2);
    for a=1:ntrials
        curr_trial=squeeze(pca_trials_shock(:,a,:));
        pca_trials_shock_omit=pca_trials_shock;
        pca_trials_shock_omit(:,a,:)=[];
        mpca_shock_omit=squeeze(mean(pca_trials_shock_omit,2));
        d_opp=sqrt( (mpca_shock_omit(1,:)-curr_trial(1,:)).^2 + (mpca_shock_omit(2,:)-curr_trial(2,:)).^2 + (mpca_shock_omit(3,:)-curr_trial(3,:)).^2 );
        d_same=sqrt( (mpca_reward(1,:)-curr_trial(1,:)).^2 + (mpca_reward(2,:)-curr_trial(2,:)).^2 + (mpca_reward(3,:)-curr_trial(3,:)).^2 );
        curr_shock=[curr_shock; (d_opp-d_same)./(d_opp+d_same)];
    end
    shock_si(z,:)=mean(curr_shock);
end
all_si=[reward_si; shock_si];

m_reward_si=mean(reward_si);
m_shock_si=mean(shock_si);
m_all_si=mean(all_si);
sem_all_si=std(all_si)./sqrt(size(all_si,1));
sem_reward_si=std(reward_si)./sqrt(size(reward_si,1));
sem_shock_si=std(shock_si)./sqrt(size(shock_si,1));
figure;
plot(time_trial,reward_si,'b');
hold on
plot(time_trial,shock_si,'r');
shadedErrorBar(time_trial,m_reward_si,sem_reward_si,'b');
hold on
shadedErrorBar(time_trial,m_shock_si,sem_shock_si,'r');
xlim([-1 5]);
ylim([-1 1]);
box off;
plot([1 1],[-1 1],'--k');
plot([0 0],[-1 1],'--k');

% figure;
% plot(time_trial,all_si,'color',[0.7 0.7 0.7]);
% hold on
% shadedErrorBar(time_trial,m_all_si,sem_all_si,'k');
% xlim([-1 5]);
% ylim([-1 1]);
% box off;
% plot([1 1],[-1 1],'--k');
% plot([0 0],[-1 1],'--k');

figure;
shadedErrorBar(time_trial,mean(euc_dist),std(euc_dist)/sqrt(size(euc_dist,1)),'k')
hold on;
box off;
xlim([-1 5])
ylim([0.005 0.02])
plot([0 0],[-0.02 0.02],'--k');
plot([1 1],[-0.02 0.02],'--k');


% m_pcs=mean(pcs_dist,2);
% sem_pcs=std(pcs_dist,0,2)/sqrt(size(pcs_dist,2));
% figure;
% subplot(3,1,1);
% shadedErrorBar(time_trial,m_pcs(1,:),sem_pcs(1,:),'k');
% hold on
% xlim([-1 5]);
% ylim([-0.01 0.015])
% plot([0 0],[-0.02 0.02],'--k');
% plot([1 1],[-0.02 0.02],'--k');
% box off
% subplot(3,1,2);
% shadedErrorBar(time_trial,m_pcs(2,:),sem_pcs(2,:),'k');
% hold on
% xlim([-1 5]);
% ylim([-0.01 0.015])
% plot([0 0],[-0.02 0.02],'--k');
% plot([1 1],[-0.02 0.02],'--k');
% box off
% subplot(3,1,3);
% shadedErrorBar(time_trial,m_pcs(3,:),sem_pcs(3,:),'k');
% hold on
% xlim([-1 5]);
% ylim([-0.01 0.015])
% plot([0 0],[-0.02 0.02],'--k');
% plot([1 1],[-0.02 0.02],'--k');
% box off
