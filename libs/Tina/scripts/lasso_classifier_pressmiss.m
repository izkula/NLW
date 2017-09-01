function [cells_lasso,weights_lasso]=lasso_classifier_pressmiss(pathname,filedate,filename)

load(strcat(pathname,filedate,'/results/',filename,'.mat'));

% for a=1:size(dfof_trials_pressed_tone,1)
%     for b=1:size(dfof_trials_pressed_tone,2)
%         dfof_trials_pressed_tone(a,b,:)=smooth(dfof_trials_pressed_tone(a,b,:),4);
%         %     maxval=max(dfof_trials_pressed_tone(a,:));
%         %     dfof_trials_pressed_tone(a,:)=dfof_trials_pressed_tone(a,:)./maxval;
%     end
% end
% 
% for a=1:size(dfof_trials_missed_tone,1)
%     for b=1:size(dfof_trials_missed_tone,2)
%         dfof_trials_missed_tone(a,b,:)=smooth(dfof_trials_missed_tone(a,b,:),4);
%         %     maxval=max(dfof_trials_missed_tone(a,:));
%         %     dfof_trials_missed_tone(a,:)=dfof_trials_missed_tone(a,:)./maxval;
%     end
% end


num_cells=size(dfof_trials_pressed_tone,1);
num_pressedtrials=size(dfof_trials_pressed_tone,2);
num_missedtrials=size(dfof_trials_missed_tone,2);
numtrials=num_pressedtrials+num_missedtrials;

t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial<=5,1,'last');

% create "X" variable for regression, only for t=0 to 5 of tone.
pressed_og=dfof_trials_pressed_tone(:,:,t_0:t_5);
pressed=reshape(pressed_og,[num_cells,size(pressed_og,2)*size(pressed_og,3)]);
missed_og=dfof_trials_missed_tone(:,:,t_0:t_5);
missed=reshape(missed_og,[num_cells,size(missed_og,2)*size(missed_og,3)]);

% % create "X" variable for regression, only for t=0 to 5 of tone.
% pressed=mean(dfof_trials_pressed_tone(:,:,t_0:t_5),3);
% missed=mean(dfof_trials_missed_tone(:,:,t_0:t_5),3);

trial_type=[zeros(1,size(pressed,2)) ones(1,size(missed,2))];
all_trials=[pressed';missed'];

cp=cvpartition(trial_type,'k',10); 

display('Calculating lasso GLM...');
[B,FitInfo]=lassoglm(all_trials,trial_type','binomial','NumLambda',25,'CV',cp,'DFMax',25);
lassoPlot(B,FitInfo,'PlotType','CV');
lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');
indx=FitInfo.Index1SE;
weights_lasso=B(:,indx);
cells_lasso=find(weights_lasso~=0);
weights_lasso=weights_lasso(cells_lasso);

savename=strcat(pathname,filedate,'/results/',filename,'_lasso_pressmiss.mat');
save(savename,'weights_lasso','cells_lasso','B','FitInfo');
end

