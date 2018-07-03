clear; clc; close all
init_EmilyS

cd '/home/svesuna/2pdata/BpodData/'

d = dir('**/*.mat')

for i = 1:numel(d)
    if regexp(d(i).name(1:4),'m830') 
try
        d_split = split(d(i).name,'_')
        name.mouse = d_split{1};;
        name.exp = d_split{2};
        name.date = [d_split{3} '_' d_split{4}];
        name.sess = d_split{5};
        
        
        %load 2p info things
        load(fullfile(d(i).folder, d(i).name))
    
        %load bpod data
       filename = dir(fullfile('/home/svesuna/Dropbox/Bpod/Data',...
           name.mouse ,name.exp,'Session Data'));
       
       %SessionData
       load(fullfile(filename(3).folder, filename(3).name))
        
       
       lickFrames = []
       figure; hold on;
       for j = 1:SessionData.nTrials
           try
               hold on;
           scatter(SessionData.RawEvents.Trial{j}.Events.Port1In,...
             j*ones(length(SessionData.RawEvents.Trial{j}.Events.Port1In),1),'k.')        
      
        %convert licks to frame numbers
         lickFrames = [lickFrames; round(...
             frameSampleRate*SessionData.RawEvents.Trial{j}.Events.Port1In + info.frame(j))];
           end
       end
     
       
       %clear unncessary things
       clearvars -except d frameSampleRate SessionData info lickFrames name  
       
      
       
       
       load(fullfile('/home/svesuna/2pdata/BpodData/',...
           name.mouse,name.exp,name.date,name.sess,...
           [name.mouse '_' name.exp '_' name.date '_' name.sess],'neuron.mat'))
       
       
       save(fullfile('/home/svesuna/Dropbox/2p_Habenula_Shared/',name.mouse, [name.mouse '_' name.exp '_' name.date '.mat']))
       
end
    end
    
end
