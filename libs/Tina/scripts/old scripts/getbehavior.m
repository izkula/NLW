close all
clear all

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/20150525/';
fn_base='thcre_trainingday3';

% open dfof file
load(strcat(pathname,'results/',fn_base,'_dfof.mat'));

% open nidaq file
fid2=fopen(strcat(pathname,'binfiles/',fn_base,'.bin'),'r');
[data,count]=fread(fid2,[6,inf],'double');
fclose(fid2);

% assign data 
timeframes=data(1,:);
camframes=data(2,:);
toneframes=data(3,:);
rewardframes=data(5,:);
lickframes=data(4,:);
shockframes=data(6,:);

% find the start inds of each 2p frame
diff_cam=diff(camframes);
start_inds=find(diff_cam<-2);
check=diff(start_inds);
del=find(check<1500)+1;
start_inds(del)=[];

% calculate time from start to end of recording for 2p frames
time=timeframes(start_inds);
time=time-time(1);

% calculate time from start to end of recording for stimulus
time_stimulus=timeframes(start_inds(1):start_inds(end));
time_stimulus=time_stimulus-time_stimulus(1);

% calculate stimulus
tone=toneframes(start_inds(1):start_inds(end));
reward=rewardframes(start_inds(1):start_inds(end));
lick=lickframes(start_inds(1):start_inds(end));
shock=shockframes(start_inds(1):start_inds(end));

% calculate starts of tones
tone_start=find(diff(tone)>2);
check=diff(tone_start);
del=find(check<2800)+1;

% finds starts of pulsed tone
tone_shock_start=tone_start(del([1:3:length(del)])-1);
tone_reward_start=tone_start;
tone_reward_start([del del([1:3:length(del)])-1])=[];

tone_shock_start(end)=[];
%tone_shock_start(1)=[];
tone_reward_start(end)=[];

cells_reward=zeros(size(dfof,1),length(tone_reward_start),42);
for a=1:length(tone_reward_start)
    
    startind=find(time>=time_stimulus(tone_reward_start(a)-20001),1,'first');
    endind=startind+41;
    
    cells_reward(:,a,:)=dfof(:,startind:endind);
end
time_reward=time(startind:endind);
time_reward=time_reward-time_reward(1)-2;

cells_shock=zeros(size(dfof,1),length(tone_shock_start),42);
for a=1:length(tone_shock_start)
    
    startind=find(time>=time_stimulus(tone_shock_start(a)-20001),1,'first');
    endind=startind+41;
    
    cells_shock(:,a,:)=dfof(:,startind:endind);
end
time_shock=time(startind:endind);
time_shock=time_shock-time_shock(1)-2;

figure(1);
for a=1:size(dfof,1)
    cell_color=rand(1,3);
    
    subplot(2,2,1)
    plot(time_reward,mean(squeeze(cells_reward(a,:,:))),'color',cell_color);%+(a-1)*0.25);
    title('Reward')
    
    
    subplot(2,2,2)
    plot(time_reward,mean(squeeze(cells_shock(a,:,:))),'color',cell_color)%+(a-1)*0.25);
    title('Shock')
    
    subplot(2,2,3)
    title('Reward')
    for b=1:length(tone_reward_start)
        plot(time_reward,squeeze(cells_reward(a,b,:)),'color',cell_color);
        hold on
    end
    hold off
    
    subplot(2,2,4)
    title('Shock')
    for b=1:length(tone_shock_start)
        plot(time_reward,squeeze(cells_shock(a,b,:)),'color',cell_color);
        hold on
    end
    hold off
    
    pause;
end

save_name=strcat(pathname,'results/',fn_base,'_stimulus.mat');
save(save_name,'time','time_stimulus','tone_reward_start','tone_shock_start');