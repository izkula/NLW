close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

class_acc=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    load(strcat(pathname,date,'/results/',filename,'_pca.mat'));
    PCA_generaterand_classifier_shuff(pathname,date,filename)
    
    num_rewardtrials=size(pca_trials_reward,2);
    
    for b=1:num_rewardtrials
        curr_trial=squeeze(pca_trials_reward(:,b,:));
        pca_trials_reward_omit=pca_trials_reward;
        pca_trials_reward_omit(:,b,:)=[];
        mpca_reward_omit=squeeze(mean(pca_trials_reward_omit,2));
        
        d_opp=sqrt( (mpca_shock(1,:)-curr_trial(1,:)).^2 + (mpca_shock(2,:)-curr_trial(2,:)).^2 + (mpca_shock(3,:)-curr_trial(3,:)).^2 );
        d_same=sqrt( (mpca_reward_omit(1,:)-curr_trial(1,:)).^2 + (mpca_reward_omit(2,:)-curr_trial(2,:)).^2 + (mpca_reward_omit(3,:)-curr_trial(3,:)).^2 );
    
        temp=zeros(1,size(pca_trials_reward,3));
        for c=1:size(pca_trials_reward,3)
            if d_same(c)<d_opp(c)
                temp(c)=1;
            else
                temp(c)=0;
            end
        end
        class_acc=[class_acc; temp];
    end
    
    
    num_shocktrials=size(pca_trials_shock,2);
    
    for b=1:num_shocktrials
        curr_trial=squeeze(pca_trials_shock(:,b,:));
        pca_trials_shock_omit=pca_trials_shock;
        pca_trials_shock_omit(:,b,:)=[];
        mpca_shock_omit=squeeze(mean(pca_trials_shock_omit,2));
        
        d_same=sqrt( (mpca_shock_omit(1,:)-curr_trial(1,:)).^2 + (mpca_shock_omit(2,:)-curr_trial(2,:)).^2 + (mpca_shock_omit(3,:)-curr_trial(3,:)).^2 );
        d_opp=sqrt( (mpca_reward(1,:)-curr_trial(1,:)).^2 + (mpca_reward(2,:)-curr_trial(2,:)).^2 + (mpca_reward(3,:)-curr_trial(3,:)).^2 );
    
        temp=zeros(1,size(pca_trials_shock,3));
        for c=1:size(pca_trials_shock,3)
            if d_same(c)<d_opp(c)
                temp(c)=1;
            else
                temp(c)=0;
            end
        end
        class_acc=[class_acc; temp];
    end

end

figure;
shadedErrorBar(time_trial,mean(class_acc),std(class_acc)./sqrt(size(class_acc,1)));
hold on
xlim([0 5]);
ylim([0 1]);
plot([0 5],[0.5 0.5],'--k');
plot([1 1],[0 1],'--k');


%%
% shuffles
numshuff=1000;
class_shuff=zeros(numshuff,size(pca_trials_reward,3));
for p=1:numshuff
    class_acc_shuff=[];
    for z=1:length(dates)
        date=dates{z};
        filename=filenames{z};

        load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
        load(strcat(pathname,date,'/results/',filename,'_pca_forshuff.mat'));
        
        % randomly assign reward and shock trials
        numtrials=size(pca_trials_pressed,2);
        shock_trials=sort(randperm(numtrials,numshock));
        reward_trials=[1:numtrials];
        reward_trials=reward_trials(~ismember(reward_trials,shock_trials));
        
        % get the reward trials closest to the shock trials
        search_reward=reward_trials;
        reward_nearshock=[];
        for a=1:length(shock_trials)
            if shock_trials(a)<search_reward(1)
                reward_nearshock(a)=search_reward(find(search_reward>shock_trials(a),1,'first'));
                search_reward(find(search_reward>shock_trials(a),1,'first'))=[];
            else
                reward_nearshock(a)=search_reward(find(search_reward<shock_trials(a),1,'last'));
                search_reward(find(search_reward<shock_trials(a),1,'last'))=[];
            end
        end

        pca_trials_reward=pca_trials_pressed(:,reward_nearshock,:);
        mpca_reward=squeeze(mean(pca_trials_reward,2));
        pca_trials_shock=pca_trials_pressed(:,shock_trials,:);
        mpca_shock=squeeze(mean(pca_trials_shock,2));
        
        num_rewardtrials=size(pca_trials_reward,2);
        for b=1:num_rewardtrials
            curr_trial=squeeze(pca_trials_reward(:,b,:));
            pca_trials_reward_omit=pca_trials_reward;
            pca_trials_reward_omit(:,b,:)=[];
            mpca_reward_omit=squeeze(mean(pca_trials_reward_omit,2));

            d_opp=sqrt( (mpca_shock(1,:)-curr_trial(1,:)).^2 + (mpca_shock(2,:)-curr_trial(2,:)).^2 + (mpca_shock(3,:)-curr_trial(3,:)).^2 );
            d_same=sqrt( (mpca_reward_omit(1,:)-curr_trial(1,:)).^2 + (mpca_reward_omit(2,:)-curr_trial(2,:)).^2 + (mpca_reward_omit(3,:)-curr_trial(3,:)).^2 );

            temp=zeros(1,size(pca_trials_reward,3));
            for c=1:size(pca_trials_reward,3)
                if d_same(c)<d_opp(c)
                    temp(c)=1;
                else
                    temp(c)=0;
                end
            end
            class_acc_shuff=[class_acc_shuff; temp];
        end


        num_shocktrials=size(pca_trials_shock,2);
        for b=1:num_shocktrials
            curr_trial=squeeze(pca_trials_shock(:,b,:));
            pca_trials_shock_omit=pca_trials_shock;
            pca_trials_shock_omit(:,b,:)=[];
            mpca_shock_omit=squeeze(mean(pca_trials_shock_omit,2));

            d_same=sqrt( (mpca_shock_omit(1,:)-curr_trial(1,:)).^2 + (mpca_shock_omit(2,:)-curr_trial(2,:)).^2 + (mpca_shock_omit(3,:)-curr_trial(3,:)).^2 );
            d_opp=sqrt( (mpca_reward(1,:)-curr_trial(1,:)).^2 + (mpca_reward(2,:)-curr_trial(2,:)).^2 + (mpca_reward(3,:)-curr_trial(3,:)).^2 );

            temp=zeros(1,size(pca_trials_shock,3));
            for c=1:size(pca_trials_shock,3)
                if d_same(c)<d_opp(c)
                    temp(c)=1;
                else
                    temp(c)=0;
                end
            end
            class_acc_shuff=[class_acc_shuff; temp];
        end

    end
    class_shuff(p,:)=mean(class_acc_shuff);
    p
end


classifier_accuracy=mean(class_acc);
pval=[];
for z=1:size(class_shuff,2)
   pval(z)=length(find(class_shuff(:,z)>=classifier_accuracy(z)))/numshuff; 
end

% shadedErrorBar(time_trial,mean(class_shuff),std(class_shuff)./sqrt(size(class_shuff,1)));

figure; plot(time_trial,pval)
hold on
plot([time_trial(1) time_trial(end)],[0.05 0.05],'--k')
xlim([0 5]);
plot([1 1],[0 1],'--k');

