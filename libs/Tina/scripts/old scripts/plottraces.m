%% PLOT SHIT

close all

load(strcat(pathname,'results/',fn_base,'_stimulus.mat'));
load(strcat(pathname,'results/',fn_base,'_dfof.mat'));

stimulus(stimulus<0)=0;

figure;
plot(time_stimulus,stimulus*(size(dfof,1)/max(stimulus)),'r');
hold on
for a=1:size(dfof,1)
    colorcode=rand(1,3);
    plot(time,dfof(a,:)+(a-1),'color',colorcode,'linewidth',1);
    hold on
end

xlabel('Time (s)');
ylabel('dF/F');
