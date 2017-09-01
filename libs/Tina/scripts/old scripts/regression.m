close all
clear all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
filedate='20160730';
fn='GRIN_pfcvta1_tonerewardshock_day4';

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/results/',fn,'.mat'));
load(strcat(currpath,'/behavior/',fn,'.mat'));

trial_type=dfof_trials_shock_press;

numtrials=size(trial_type,2);
numtimepoints=length(time_trial);
numcells=size(trial_type,1);

t_1=find(time_trial>=1,1,'first');
t_2=find(time_trial>=2,1,'first');
behav=zeros(1,length(time_trial));
behav(t_1:t_2)=5;

% convolve behavior regressor with exponential
tau=0.25;
cirf=exp(-time_trial/tau);
cirf=cirf/sum(cirf);
behav_conv=conv(behav,cirf);
behav_conv=behav_conv(1:length(cirf));
behav_regress=repmat(behav_conv,1,numtrials);

cells_type=[];
for a=1:numcells
    cells_type(a,:)=reshape(squeeze(trial_type(a,:,:))',1,numtimepoints*numtrials);
    mcells_type(a,:)=mean(squeeze(trial_type(a,:,:)));
    semcells_type(a,:)=std(squeeze(trial_type(a,:,:)))/sqrt(numtrials);
end

[rho,p]=corr(behav_regress',cells_type');

behav_cells=find(p<0.05);
pos=find(rho(behav_cells)>=0.2);
behav_cells=behav_cells(pos);

% for a=1:length(behav_cells)
%     figure(1);
%     plot(cells_type(behav_cells(a),:));
%     hold on;
%     plot([1:numtrials*numtimepoints],behav_regress,'r','linewidth',2);
%     hold off
%     pause;
%     
% end
%plot([1:numtrials*numtimepoints],behav_regress,'r','linewidth',2);
%plot(time_trial,behav);
display(length(behav_cells));
display(length(behav_cells)/numcells);

trial_type2=dfof_trials_reward_press;
numtrials2=size(trial_type2,2);
cells_type2=[];
for a=1:numcells
    cells_type2(a,:)=reshape(squeeze(trial_type2(a,:,:))',1,numtimepoints*numtrials2);
    mcells_type2(a,:)=mean(squeeze(trial_type2(a,:,:)));
    semcells_type2(a,:)=std(squeeze(trial_type2(a,:,:)))/sqrt(numtrials2);
end
behav_regress2=repmat(behav_conv,1,numtrials2);

[rho2,p2]=corr(behav_regress2',cells_type2');
% overlap=find(p2(behav_cells)<0.05);


t_0=find(time_trial>=0,1,'first');
t_5=find(time_trial>=5,1,'first');
behav2=zeros(1,length(time_trial));
behav2(t_0:t_5)=5;

trial_type3=dfof_trials_missed_tone;
numtrials3=size(trial_type3,2);
cells_type3=[];
for a=1:numcells
    cells_type3(a,:)=reshape(squeeze(trial_type3(a,:,:))',1,numtimepoints*numtrials3);
    mcells_type3(a,:)=mean(squeeze(trial_type3(a,:,:)));
    semcells_type3(a,:)=std(squeeze(trial_type3(a,:,:)))/sqrt(numtrials3);
end


% only look at trials where he pressed within 1s
pressed_trials=sort([reward_trials shock_trials]);
reaction_times=time_behav(lever_retract(pressed_trials))-time_behav(lever_extend(pressed_trials));
figure; plot([1:length(reaction_times)],reaction_times,'o');
del_long=find(reaction_times>1);
hold on
plot(del_long,reaction_times(del_long),'x');
shock_inds=ismember(pressed_trials,shock_trials);
shock_inds=find(shock_inds==1);
plot(shock_inds,reaction_times(shock_inds),'og');

trial_type4=dfof_trials_pressed_tone;
trial_type4(:,del_long,:)=[];
numtrials4=size(trial_type4,2);
cells_type4=[];
for a=1:numcells
    cells_type4(a,:)=reshape(squeeze(trial_type4(a,:,:))',1,numtimepoints*numtrials4);
    mcells_type4(a,:)=mean(squeeze(trial_type4(a,:,:)));
    semcells_type4(a,:)=std(squeeze(trial_type4(a,:,:)))/sqrt(numtrials4);
end



% for a=1:length(behav_cells)
%     figure(3);
%     %shadedErrorBar(time_trial,mcells_type(behav_cells(a),:),semcells_type(behav_cells(a),:));
%     plot(time_trial,mcells_type(behav_cells(a),:),'r');
%     hold on;
%     plot(time_trial,mcells_type2(behav_cells(a),:),'b');
%     plot(time_trial,behav,'k');
%     xlim([-2 5]);
%     hold off
%     
%     figure(4);
%     %shadedErrorBar(time_trial,mcells_type(behav_cells(a),:),semcells_type(behav_cells(a),:));
%     plot(time_trial,mcells_type3(behav_cells(a),:),'r');
%     hold on;
%     plot(time_trial,mcells_type4(behav_cells(a),:),'b');
%     plot(time_trial,behav2,'k');
%     xlim([-2 5]);
%     hold off
%     
%     pause;
% end

trial_type5=dfof_trials_pressed_tone;
peaks=max(trial_type5(:,:,t_0:t_5),[],3);

for a=1:length(behav_cells)
    figure(3);
%     for b=1:numtrials
%         plot(time_trial,squeeze(trial_type(behav_cells(a),b,:)),'r');
%         hold on;
%     end
%     for b=1:numtrials2
%         plot(time_trial,squeeze(trial_type2(behav_cells(a),b,:)),'b');
%     end
    shadedErrorBar(time_trial,mcells_type(behav_cells(a),:),semcells_type(behav_cells(a),:),'r');
    hold on
    shadedErrorBar(time_trial,mcells_type2(behav_cells(a),:),semcells_type2(behav_cells(a),:),'b');
    plot(time_trial,behav,'k');
    xlim([-2 5]);
    title(strcat('Cell: ',num2str(behav_cells(a)),'Rho shock:',num2str(rho(behav_cells(a))),'Rho reward:',num2str(rho2(behav_cells(a))),'P reward:',num2str(p2(behav_cells(a)))));
    hold off
    
    figure(4);
    for b=1:numtrials3
        plot(time_trial,smooth(squeeze(trial_type3(behav_cells(a),b,:)),5),'r');
        hold on;
    end
    for b=1:numtrials4
        plot(time_trial,smooth(squeeze(trial_type4(behav_cells(a),b,:)),5),'b');
    end
    plot(time_trial,behav2,'k');
    xlim([-2 5]);
        title(strcat('Rho shock:',num2str(rho(behav_cells(a))),'Rho reward:',num2str(rho2(behav_cells(a))),'P reward:',num2str(p2(behav_cells(a)))));

    hold off
    
    figure(5);
    plot(reaction_times,peaks(behav_cells(a),:),'o');
    hold off
    
    pause;
end
