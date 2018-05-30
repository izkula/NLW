clear; clc; close all
init_samz
%
f = dir(fullfile('~/2pdata/20180417/','*'));

%%
 j = 1

for i = 1:numel(f)
    try
        %get concatenated file names
        if contains(f(i).name,'fp')
            filepaths{j} = fullfile(f(i).folder, f(i).name);
            folderpaths{j} = f(i).folder;

            j = j+1;
        end
    end
end
%%

 %Define and reset important variables of metadata

datasets = {}; 
for ii = 1:numel(filepaths)
s = split(filepaths{ii},'/')
name = split(s{end},'_')
name = name{1}

%load NIDAQ information
nidaqSampelRate = 5000;
try
[t,d] = LoadNidaqOutput(fullfile(filepaths{ii}, 'nidaq.mat'),5);

%channels
%1  - stim
%2 -trial start end
%3 - 407
%4 - 405
%5 - running
figure;
for k = 1:5
    subplot(5,1,k)
    plot(d(k,:))
end


%EXTRACT TRIAL START AND ENDS
meta.raw.sync = d(2,:) %save raw trial trace
trialTicks = find(d(2,:)> max(d(2,:)*0.98)); 

%^SAVE OUT FOR TRIAL START AND ENDS
meta.trialStart = trialTicks(1);
meta.trialEnd = trialTicks(end);

%EXTRACT RUN SPEEDS
runTicks = d(5,:) < -0.1;
tickBin = nidaqSampelRate / 20; %for a bin of 500 ms

runCounts = sum(reshape(runTicks, tickBin, []));
runCountsSmooth = smooth(runCounts, 50); 

%SAVE OUT FPR RUNNING
meta.raw.run = d(2,:)
meta.runCounts = runCounts;
meta.runCountsSmooth = runCountsSmooth;
meta.runTickBin = tickBin
meta.nidaqSampleRate = nidaqSampelRate

%EXTRACT FP SIGNAL
Y = smooth(d(3,:),1000);
XX = smooth(d(4,:),1000);
bls = polyfit(XX,Y,1);
Y_fit = bls(1) * XX + bls(2);
Y_dF = Y - Y_fit;

%SAVE OUT FP
meta.fp = Y_dF
meta.raw.fpsig = d(3,:);
meta.raw.fpcon = d(4,:);



if~exist(fullfile('~', '2presults', '20180417', name))
    mkdir(fullfile('~', '2presults', '20180417', name))
end


save(fullfile('~', '2presults', '20180417', name, 'metadata.mat'),'meta')

end

    
end
