function run_loadbinfile(pathname,filedate,fn,frametime)

display('Opening csv file...');

currpath=strcat(pathname,filedate);


% open binfile
fid=fopen(strcat(currpath,'/binfiles/',fn,'.csv'),'r');
M=textscan(fid,'%f32, %f32, %f32, %f32, %f32, %f32, %f32','HeaderLines',1);
%close(fid);
%M{7}(end+1)=0;
M=cat(2,M{:});
M=M';

samplerate=50000;

% assign data
time_behav=M(1,:)/10e2;

% open dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof');
%time=[0:frametime:(length(dfof)-1)*frametime];

frames=M(2,:);
frame_start=find(diff(frames)>1.8)+1;
del=find(diff(frame_start)<0.1*samplerate)+1;
frame_start(del)=[];

figure; plot(diff(frame_start),'o');
title('detected 2p frames');

time=time_behav(frame_start);

% delete extra behavior
if time_behav(end)>time(end)
    del=find(time_behav>=time(end));
    M(:,del)=[];
    time_behav(del)=[];
end

tone=M(3,:);
lever=M(4,:);
reward=M(5,:);
lick=M(6,:);
shock=M(7,:);

% calculate starts of lever
tone_start=find(diff(tone)>1.8)+1;
lever_extend=find(diff(lever)<-1.8)+1;
lever_retract=find(diff(lever)>1.8)+1;
reward_start=find(diff(reward)>1.8)+1;
shock_start=find(diff(shock)>1.8)+1;

% delete extra 
del=find(diff(reward_start)<10*samplerate)+1;
reward_start(del)=[];
del=find(diff(tone_start)<10*samplerate)+1;
tone_start(del)=[];
del=find(diff(shock_start)<10*samplerate)+1;
shock_start(del)=[];
del=find(diff(lever_extend)<10*samplerate)+1;
lever_extend(del)=[];
del=find(diff(lever_retract)<10*samplerate)+1;
lever_retract(del)=[];

% delete extra tone/extend
if length(lever_extend)>length(lever_retract)
    del=[length(lever_retract)+1:length(lever_extend)];
    lever_extend(del)=[];
    tone_start(del)=[];
end

% delete spurious noise
del=[];
allshock_end=find(diff(shock)<-1.8)+1;
for a=1:length(shock_start)
    shock_end=allshock_end(find(allshock_end>shock_start(a),1,'first'));
    if shock_end-shock_start(a)<0.9*samplerate
       del=[del a];
    end
end
shock_start(del)=[];
del=[];
allreward_end=find(diff(reward)<-1.8)+1;
for a=1:length(reward_start)
    reward_end=allreward_end(find(allreward_end>reward_start(a),1,'first'));
    if reward_end-reward_start(a)<0.24*samplerate
       del=[del a];
    end
end
reward_start(del)=[];

display(strcat('Tone start: ',num2str(length(tone_start))));
display(strcat('Lever extend: ',num2str(length(lever_extend))));
display(strcat('Lever retract: ',num2str(length(lever_retract))));
display(strcat('Reward: ',num2str(length(reward_start))));
display(strcat('Shock: ',num2str(length(shock_start))));
    
% find successful presses
reward_trials=[];
missed_trials=[];
shock_trials=[];
for a=1:length(lever_retract)
    search_reward=find(reward_start>=lever_retract(a),1,'first');
    search_shock=find(shock_start>=lever_retract(a),1,'first');
    if reward_start(search_reward)-lever_retract(a)<5*samplerate
        reward_trials=[reward_trials a];
    elseif shock_start(search_shock)-lever_retract(a)<5*samplerate
        shock_trials=[shock_trials a];
    else
        missed_trials=[missed_trials a];
    end
end

display(strcat('# reward trials: ',num2str(length(reward_trials))));
display(strcat('# shock trials: ',num2str(length(shock_trials))));
display(strcat('# missed trials: ',num2str(length(missed_trials))));
display(strcat('# total trials: ',num2str(length(missed_trials)+length(reward_trials)+length(shock_trials))));


if ~exist(strcat(currpath,'/behavior/'),'dir')
    mkdir(strcat(currpath,'/behavior/'));
end

savename=strcat(currpath,'/behavior/',fn,'.mat');
save(savename,'tone_start','lever_extend','lever_retract','reward_start','shock_start','reward_trials','missed_trials','shock_trials','time_behav','lick','time');
