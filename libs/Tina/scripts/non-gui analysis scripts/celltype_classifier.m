function celltype_classifier(pathname,filedate,filename)

load(strcat(pathname,filedate,'/results/',filename,'.mat'));

num_cells=size(dfof_trials_reward_press,1);
num_rewardtrials=size(dfof_trials_reward_press,2);
num_shocktrials=size(dfof_trials_shock_press,2);
numtrials=num_rewardtrials+num_shocktrials;

cells_class=zeros(num_cells,numtrials,length(time_trial));
for d=1:num_cells
    class_acc=[];
    for b=1:num_rewardtrials
        curr_trial=squeeze(dfof_trials_reward_press(d,b,:))';
        dfof_trials_reward_omit=squeeze(dfof_trials_reward_press(d,:,:));
        dfof_trials_reward_omit(b,:)=[];
        mtrials_reward_omit=mean(dfof_trials_reward_omit,1);
        
        d_opp=sqrt( (mtrials_shock_press(d,:)-curr_trial).^2 );
        d_same=sqrt( (mtrials_reward_omit-curr_trial).^2 );
        
        temp=zeros(1,size(dfof_trials_reward_press,3));
        for c=1:size(dfof_trials_reward_press,3)
            if d_same(c)<d_opp(c)
                temp(c)=1;
            else
                temp(c)=0;
            end
        end
        class_acc=[class_acc; temp];
    end
    
    for b=1:num_shocktrials
        curr_trial=squeeze(dfof_trials_shock_press(d,b,:))';
        dfof_trials_shock_omit=squeeze(dfof_trials_shock_press(d,:,:));
        dfof_trials_shock_omit(b,:)=[];
        mtrials_shock_omit=mean(dfof_trials_shock_omit,1);
        
        d_same=sqrt( (mtrials_shock_omit-curr_trial).^2 );
        d_opp=sqrt( (mtrials_reward_press(d,:)-curr_trial).^2 );
        
        temp=zeros(1,size(dfof_trials_shock_press,3));
        for c=1:size(dfof_trials_shock_press,3)
            if d_same(c)<d_opp(c)
                temp(c)=1;
            else
                temp(c)=0;
            end
        end
        class_acc=[class_acc; temp];
    end
    cells_class(d,:,:)=class_acc;
end


% shuffles
numshuff=1000;
cells_class_shuff=zeros(d,numshuff,length(time_trial));
for d=1:num_cells
    dfof_trials_all_press=cat(2,dfof_trials_reward_press,dfof_trials_shock_press);
    dfof_trials_all_press=squeeze(dfof_trials_all_press(d,:,:));
    
    for p=1:numshuff
        class_acc_shuff=[];
        
        % randomly assign reward and shock trials
        shock_trials_choose=sort(randperm(numtrials,num_shocktrials));
        reward_trials_choose=[1:numtrials];
        reward_trials_choose=reward_trials_choose(~ismember(reward_trials_choose,shock_trials_choose));
        
        dfof_trials_reward_shuff=dfof_trials_all_press(reward_trials_choose,:);
        mtrials_reward_shuff=mean(dfof_trials_reward_shuff,1);
        dfof_trials_shock_shuff=dfof_trials_all_press(shock_trials_choose,:);
        mtrials_shock_shuff=mean(dfof_trials_shock_shuff,1);
        
        for b=1:num_rewardtrials
            curr_trial=dfof_trials_reward_shuff(b,:);
            dfof_trials_reward_omit_shuff=dfof_trials_reward_shuff;
            dfof_trials_reward_omit_shuff(b,:)=[];
            mtrials_reward_omit_shuff=mean(dfof_trials_reward_omit_shuff,1);
            
            d_opp=sqrt( (mtrials_shock_shuff-curr_trial).^2 );
            d_same=sqrt( (mtrials_reward_omit_shuff-curr_trial).^2 );
            
            temp=zeros(1,length(time_trial));
            for c=1:length(time_trial)
                if d_same(c)<d_opp(c)
                    temp(c)=1;
                else
                    temp(c)=0;
                end
            end
            class_acc_shuff=[class_acc_shuff; temp];
        end
        
        for b=1:num_shocktrials
            curr_trial=dfof_trials_shock_shuff(b,:);
            dfof_trials_shock_omit_shuff=dfof_trials_shock_shuff;
            dfof_trials_shock_omit_shuff(b,:)=[];
            mtrials_shock_omit_shuff=mean(dfof_trials_shock_omit_shuff,1);
            
            d_same=sqrt( (mtrials_shock_omit_shuff-curr_trial).^2 );
            d_opp=sqrt( (mtrials_reward_shuff-curr_trial).^2 );
            
            temp=zeros(1,length(time_trial));
            for c=1:length(time_trial)
                if d_same(c)<d_opp(c)
                    temp(c)=1;
                else
                    temp(c)=0;
                end
            end
            class_acc_shuff=[class_acc_shuff; temp];
        end
        
        cells_class_shuff(d,p,:)=mean(class_acc_shuff);
    end
    d
end

t_1=find(time_trial>=1,1,'first');
t_3=find(time_trial>=3,1,'first');
m_cells_class=squeeze(mean(cells_class,2));
classifier_accuracy=mean(m_cells_class(:,t_1:t_3),2);
pval=[];
for z=1:size(cells_class_shuff,1)
    currcell=squeeze(cells_class_shuff(z,:,:));
    currcell=mean(currcell(:,t_1:t_3),2);
    pval(z)=length(find(currcell>=classifier_accuracy(z)))/numshuff;
end

discriminators=find(pval<0.05);

% shadedErrorBar(time_trial,mean(class_shuff),std(class_shuff)./sqrt(size(class_shuff,1)));

% for a=1:length(discriminators)
%     figure(1);
%     shadedErrorBar(time_trial,mtrials_reward_press(discriminators(a),:),semtrials_reward_press(discriminators(a),:),'b');
%     hold on
%     shadedErrorBar(time_trial,mtrials_shock_press(discriminators(a),:),semtrials_shock_press(discriminators(a),:),'r');
%     hold off;
%     pause;
% end

num_discriminators=length(discriminators);
display(num_discriminators);

savename=strcat(pathname,filedate,'/results/',filename,'_celltypes.mat');
save(savename,'discriminators','cells_class','cells_class_shuff');
end

