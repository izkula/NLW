shock_start=find(time_trial>=1,1,'first');
shock_end=find(time_trial>=4,1,'first');

base_inds=[1:shock_start-1];
shock_inds=[shock_start:shock_end];


for z=1:size(dfof_trials_shock_press,1)
    curr_cell=squeeze(dfof_trials_shock_press(z,:,:));
    
    baseline(z,:)=curr_cell(:,base_inds);
    thresh(z)=std(baseline)*2.5;
    
    for a=1:size(dfof_trials_shock_press,2)
        active_trials=find()
        
    end
    
    
    shock(z,:)=mean(curr_cell(:,shock_inds),2);
    p(z)=signrank(baseline(z,:),shock(z,:));
end

shock_cells=find(p<0.05);

for b=1:length(shock_cells)
    a=shock_cells(b);
    figure(1);
    subplot(2,4,1);
    shadedErrorBar(time_trial,mtrials_reward_press(a,:),semtrials_reward_press(a,:));
    hold on
    plot([0 0],[-0.2 0.2],'--k');
    plot([1 1],[-0.2 0.2],'--k');
    hold off
    ylim([-0.2 0.2])
    xlim([-2 10])
    title(strcat('Cell ',num2str(a),'Reward'));
    subplot(2,4,2);
    shadedErrorBar(time_trial,mtrials_shock_press(a,:),semtrials_shock_press(a,:));
    hold on
    plot([0 0],[-0.2 0.2],'--k');
    plot([1 1],[-0.2 0.2],'--k');
    plot([2 2],[-0.2 0.2],'--k');
    hold off
    ylim([-0.2 0.2])
    xlim([-2 10])
    title('Shock');
    subplot(2,4,3);
    shadedErrorBar(time_trial,mtrials_pressed_tone(a,:),semtrials_pressed_tone(a,:));
    hold on
    plot([0 0],[-0.2 0.2],'--k');
    plot([5 5],[-0.2 0.2],'--k');
    hold off
    ylim([-0.2 0.2])
    xlim([-2 5])
    title('Tone pressed');
    subplot(2,4,4);
    shadedErrorBar(time_trial,mtrials_missed_tone(a,:),semtrials_missed_tone(a,:));
    hold on
    plot([0 0],[-0.2 0.2],'--k');
    plot([5 5],[-0.2 0.2],'--k');
    hold off
    ylim([-0.2 0.2])
    xlim([-2 10])
    title('Tone missed');
    
    subplot(2,4,5);
    for b=1:size(dfof_trials_reward_press,2)
        plot(time_trial,squeeze(dfof_trials_reward_press(a,b,:)));
        ylim([-0.5 1])
        xlim([-2 10])
        hold on
    end
    plot([0 0],[-0.5 1],'--k');
    plot([1 1],[-0.5 1],'--k');
    hold off
    subplot(2,4,6);
    for b=1:size(dfof_trials_shock_press,2)
        plot(time_trial,squeeze(dfof_trials_shock_press(a,b,:)));
        ylim([-0.5 1])
        xlim([-2 10])
        hold on
    end
    plot([0 0],[-0.5 1],'--k');
    plot([1 1],[-0.5 1],'--k');
    plot([2 2],[-0.5 1],'--k');
    hold off
    subplot(2,4,7);
    for b=1:size(dfof_trials_pressed_tone,2)
        plot(time_trial,squeeze(dfof_trials_pressed_tone(a,b,:)));
        ylim([-0.5 1])
        xlim([-2 5])
        hold on
    end
    plot([0 0],[-0.5 1],'--k');
    plot([5 5],[-0.5 1],'--k');
    hold off
    subplot(2,4,8);
    for b=1:size(dfof_trials_missed_tone,2)
        plot(time_trial,squeeze(dfof_trials_missed_tone(a,b,:)));
        ylim([-0.5 1])
        xlim([-2 10])
        hold on
    end
    plot([0 0],[-0.5 1],'--k');
    plot([5 5],[-0.5 1],'--k');
    hold off
    pause;
end