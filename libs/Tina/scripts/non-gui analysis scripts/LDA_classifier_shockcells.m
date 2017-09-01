close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

% dates={'20160529','20160529','20160730','20160730','20160731'};
% filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

class_acc=[];
for z=1:length(dates)
    w_cells=[];
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/results/',filename,'.mat'));
%     load(strcat(pathname,date,'/regression/',filename,'_leveronly.mat'));
%     load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
%     load(strcat(pathname,date,'/regression/',filename,'_rewardonly.mat'));
    [cells_lasso,weights_lasso]=lasso_classifier_pressmiss_shock(pathname,date,filename);

    num_pressedtrials=size(dfof_trials_pressed_tone,2);
    num_missedtrials=size(dfof_trials_missed_tone,2);
    numtrials=num_pressedtrials+num_missedtrials;

    t_0=find(time_trial>=0,1,'first');
    t_5=find(time_trial<=5,1,'last');
 
    
    % USE ONLY THE SHOCK LASSO NEURONS
    num_cells=length(cells_lasso);
    dfof_trials_pressed_tone=dfof_trials_pressed_tone(cells_lasso,:,:);
    dfof_trials_missed_tone=dfof_trials_missed_tone(cells_lasso,:,:);

%     % USE ALL NEURONS
%     num_cells=size(dfof_trials_pressed_tone,1);
    
    % create "X" variable for regression, only for t=0 to 5 of tone.
    pressed_og=dfof_trials_pressed_tone(:,:,t_0:t_5);
    pressed=reshape(pressed_og,[num_cells,size(pressed_og,2)*size(pressed_og,3)]);
    missed_og=dfof_trials_missed_tone(:,:,t_0:t_5);
    missed=reshape(missed_og,[num_cells,size(missed_og,2)*size(missed_og,3)]);
    all_trials=[pressed'; missed'];
    trial_types=[zeros(1,size(pressed,2)) ones(1,size(missed,2))];

%     % create "X" variable for regression, only for t=0 to 5 of tone.
%     pressed=squeeze(mean(dfof_trials_pressed_tone(:,:,t_0:t_5),3));
%     missed=squeeze(mean(dfof_trials_missed_tone(:,:,t_0:t_5),3));
%     all_trials=[pressed'; missed'];
%     trial_types=[zeros(1,size(pressed,2)) ones(1,size(missed,2))];
        
    cp=cvpartition(trial_types,'k',5); 
    classf = @(XTRAIN,ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain));
    
%     class_acc(z)=1-crossval('mcr',all_trials,trial_types','predfun',classf,'leaveout',1);
    class_acc(z)=1-crossval('mcr',all_trials,trial_types','predfun',classf,'partition',cp,'mcr',10);
    Mdl = fitcdiscr(all_trials,trial_types');
    w_cells=Mdl.Coeffs(2,1).Linear;
    const=Mdl.Coeffs(2,1).Const;
    
%     savename=strcat(pathname,date,'/results/',filename,'_lda_weights.mat');
%     save(savename,'Mdl');

end

%%
% shuffles
numshuff=1000;
class_shuff=zeros(numshuff,length(dates));
for p=1:numshuff
    class_acc_shuff=[];
    for z=1:length(dates)
        date=dates{z};
        filename=filenames{z};
        
        load(strcat(pathname,date,'/results/',filename,'.mat'));
%         load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
        load(strcat(pathname,date,'/results/',filename,'_lasso_pressmiss_shock.mat'));
        
        num_pressedtrials=size(dfof_trials_pressed_tone,2);
        num_missedtrials=size(dfof_trials_missed_tone,2);
        numtrials=num_pressedtrials+num_missedtrials;
        
%         % USE ONLY THE MODULATED NEURONS
%         cells_task=sort([shock_only reward_only lever_only]);
%         num_cells=length(cells_task);
%         dfof_trials_pressed_tone=dfof_trials_pressed_tone(cells_task,:,:);
%         dfof_trials_missed_tone=dfof_trials_missed_tone(cells_task,:,:);
        
        % USE ONLY THE SHOCK LASSO NEURONS
        num_cells=length(cells_lasso);
        dfof_trials_pressed_tone=dfof_trials_pressed_tone(cells_lasso,:,:);
        dfof_trials_missed_tone=dfof_trials_missed_tone(cells_lasso,:,:);
        
        % create "X" variable for regression, only for t=0 to 5 of tone.
        pressed_og=dfof_trials_pressed_tone(:,:,t_0:t_5);
        pressed=reshape(pressed_og,[num_cells,size(pressed_og,2)*size(pressed_og,3)]);
        missed_og=dfof_trials_missed_tone(:,:,t_0:t_5);
        missed=reshape(missed_og,[num_cells,size(missed_og,2)*size(missed_og,3)]);
        all_trials=[pressed'; missed'];
        trial_types=[zeros(1,size(pressed,2)) ones(1,size(missed,2))];
        order=randperm(length(trial_types));
        trial_types=trial_types(order);

%         % create "X" variable for regression, only for t=0 to 5 of tone.
%         pressed=squeeze(mean(dfof_trials_pressed_tone(:,:,t_0:t_5),3));
%         missed=squeeze(mean(dfof_trials_missed_tone(:,:,t_0:t_5),3));
%         all_trials=[pressed'; missed'];
%         trial_types=[zeros(1,size(pressed,2)) ones(1,size(missed,2))];
%         order=randperm(length(trial_types));
%         trial_types=trial_types(order);

        cp=cvpartition(trial_types,'k',5); 
        classf = @(XTRAIN,ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain));
        
%         class_acc_shuff(z)=1-crossval('mcr',all_trials,trial_types','predfun',classf,'leaveout',1);
        class_acc_shuff(z)=1-crossval('mcr',all_trials,trial_types','predfun',classf,'partition',cp,'mcr',10);
    end
    class_shuff(p,:)=class_acc_shuff;
    p
end

for z=1:length(dates)
    pval(z)=length(find(class_shuff(:,z)>=class_acc(z)))/numshuff;
end


classifier_accuracy=mean(class_acc);
pval_all=length(find(mean(class_shuff,2)>=classifier_accuracy))/numshuff;

figure;
plot([1:5],class_acc,'.','markersize',30);
hold on
plot([0.5 5.5],[0.5 0.5],'--k');
errorbar([1:5],mean(class_shuff),1.96*std(class_shuff));
ylim([0 1])
xlim([0.5 5.5])