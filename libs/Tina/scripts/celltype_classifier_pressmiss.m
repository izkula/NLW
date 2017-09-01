function celltype_classifier_pressmiss(pathname,filedate,filename)

load(strcat(pathname,filedate,'/results/',filename,'.mat'));
load(strcat(pathname,filedate,'/behavior/',filename,'.mat'),'time_behav','lever_retract','lever_extend','shock_trials','reward_trials');

pressed_trials=sort([reward_trials shock_trials]);
reaction_times=time_behav(lever_retract(pressed_trials))-time_behav(lever_extend(pressed_trials));
del_long=find(reaction_times>1);
hold on
plot(del_long,reaction_times(del_long),'x');
shock_inds=ismember(pressed_trials,shock_trials);
shock_inds=find(shock_inds==1);
plot(shock_inds,reaction_times(shock_inds),'og');

trial_type4=dfof_trials_pressed_tone;
trial_type4(:,del_long,:)=[];
numtrials4=size(trial_type4,2);


num_cells=size(dfof_trials_pressed_tone,1);
num_pressedtrials=size(dfof_trials_pressed_tone,2);
num_missedtrials=size(dfof_trials_missed_tone,2);
numtrials=num_pressedtrials+num_missedtrials;

cells_class=zeros(num_cells,numtrials);
for d=1:num_cells
    class_acc=[];
    for b=1:num_pressedtrials
        curr_trial=squeeze(dfof_trials_pressed_tone(d,b,:))';
        dfof_trials_pressed_omit=squeeze(dfof_trials_pressed_tone(d,:,:));
        dfof_trials_pressed_omit(b,:)=[];
        mtrials_pressed_omit=mean(dfof_trials_pressed_omit,1);
        
%         d_opp=sqrt( (mtrials_missed_tone(d,:)-curr_trial).^2 );
%         d_same=sqrt( (mtrials_pressed_omit-curr_trial).^2 );
%         temp=zeros(1,size(dfof_trials_pressed_tone,3));
%         for c=1:size(dfof_trials_pressed_tone,3)
%             if d_same(c)<d_opp(c)
%                 temp(c)=1;
%             else
%                 temp(c)=0;
%             end
%         end
        
        d_opp=mean(sqrt( (mtrials_missed_tone(d,:)-curr_trial).^2 ));
        d_same=mean(sqrt( (mtrials_pressed_omit-curr_trial).^2 ));
        if d_same<d_opp
            temp=1;
        else
            temp=0;
        end
        
        class_acc=[class_acc; temp];
    end
    
    for b=1:num_missedtrials
        curr_trial=squeeze(dfof_trials_missed_tone(d,b,:))';
        dfof_trials_missed_omit=squeeze(dfof_trials_missed_tone(d,:,:));
        dfof_trials_missed_omit(b,:)=[];
        mtrials_missed_omit=mean(dfof_trials_missed_omit,1);
        
%         d_same=sqrt( (mtrials_missed_omit-curr_trial).^2 );
%         d_opp=sqrt( (mtrials_pressed_tone(d,:)-curr_trial).^2 );        
%         temp=zeros(1,size(dfof_trials_missed_tone,3));
%         for c=1:size(dfof_trials_missed_tone,3)
%             if d_same(c)<d_opp(c)
%                 temp(c)=1;
%             else
%                 temp(c)=0;
%             end
%         end

        d_same=mean(sqrt( (mtrials_missed_omit-curr_trial).^2 ));
        d_opp=mean(sqrt( (mtrials_pressed_tone(d,:)-curr_trial).^2 ));
        if d_same<d_opp
            temp=1;
        else
            temp=0;
        end

        class_acc=[class_acc; temp];
    end
    cells_class(d,:,:)=class_acc;
end


% shuffles
numshuff=1000;
cells_class_shuff=zeros(d,numshuff);
for d=1:num_cells
    dfof_trials_all_press=cat(2,dfof_trials_pressed_tone,dfof_trials_missed_tone);
    dfof_trials_all_press=squeeze(dfof_trials_all_press(d,:,:));
    
    for p=1:numshuff
        class_acc_shuff=[];
        
        % randomly assign pressed and missed trials
        missed_trials_choose=sort(randperm(numtrials,num_missedtrials));
        pressed_trials_choose=[1:numtrials];
        pressed_trials_choose=pressed_trials_choose(~ismember(pressed_trials_choose,missed_trials_choose));
        
        dfof_trials_pressed_shuff=dfof_trials_all_press(pressed_trials_choose,:);
        mtrials_pressed_shuff=mean(dfof_trials_pressed_shuff,1);
        dfof_trials_missed_shuff=dfof_trials_all_press(missed_trials_choose,:);
        mtrials_missed_shuff=mean(dfof_trials_missed_shuff,1);
        
        for b=1:num_pressedtrials
            curr_trial=dfof_trials_pressed_shuff(b,:);
            dfof_trials_pressed_omit_shuff=dfof_trials_pressed_shuff;
            dfof_trials_pressed_omit_shuff(b,:)=[];
            mtrials_pressed_omit_shuff=mean(dfof_trials_pressed_omit_shuff,1);
            
            %             d_opp=sqrt( (mtrials_missed_shuff-curr_trial).^2 );
            %             d_same=sqrt( (mtrials_pressed_omit_shuff-curr_trial).^2 );
            %             temp=zeros(1,length(time_trial));
            %             for c=1:length(time_trial)
            %                 if d_same(c)<d_opp(c)
            %                     temp(c)=1;
            %                 else
            %                     temp(c)=0;
            %                 end
            %             end
            
            d_opp=mean(sqrt( (mtrials_missed_shuff-curr_trial).^2 ));
            d_same=mean(sqrt( (mtrials_pressed_omit_shuff-curr_trial).^2 ));
            if d_same<d_opp
                temp=1;
            else
                temp=0;
            end
            
            class_acc_shuff=[class_acc_shuff; temp];
        end
        
        for b=1:num_missedtrials
            curr_trial=dfof_trials_missed_shuff(b,:);
            dfof_trials_missed_omit_shuff=dfof_trials_missed_shuff;
            dfof_trials_missed_omit_shuff(b,:)=[];
            mtrials_missed_omit_shuff=mean(dfof_trials_missed_omit_shuff,1);
            
            %             d_same=sqrt( (mtrials_missed_omit_shuff-curr_trial).^2 );
            %             d_opp=sqrt( (mtrials_pressed_shuff-curr_trial).^2 );
            %             temp=zeros(1,length(time_trial));
            %             for c=1:length(time_trial)
            %                 if d_same(c)<d_opp(c)
            %                     temp(c)=1;
            %                 else
            %                     temp(c)=0;
            %                 end
            %             end
            
            d_same=mean(sqrt( (mtrials_missed_omit_shuff-curr_trial).^2 ));
            d_opp=mean(sqrt( (mtrials_pressed_shuff-curr_trial).^2 ));
            if d_same<d_opp
                temp=1;
            else
                temp=0;
            end
            
            class_acc_shuff=[class_acc_shuff; temp];
        end
        
        cells_class_shuff(d,p)=mean(class_acc_shuff);
    end
    d
end

t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial<5,1,'last');
m_cells_class=squeeze(mean(cells_class,2));
classifier_accuracy=m_cells_class;
pval=[];
for z=1:size(cells_class_shuff,1)
    currcell=squeeze(cells_class_shuff(z,:));
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

savename=strcat(pathname,filedate,'/results/',filename,'_pressmiss.mat');
save(savename,'discriminators','cells_class','cells_class_shuff');
end

