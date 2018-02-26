direc = fullfile(''); %%% Set this if you will be using frequently

% dirname = uigetdir(['C:\2pdata\', direc])

dirnames = {
%     'C:\2pdata\BpodData\m119\ImageNLWTrainStimGoNoGo\Jul13_2017\'; 
%             'C:\2pdata\BpodData\m118';
%             'C:\2pdata\BpodData\m407';
%             'C:\2pdata\BpodData\m295';
%             'C:\2pdata\BpodData\m120';

%     'C:\2pdata\BpodData\rbpm12\ImageNLWTrainOlfGNG';
%     'C:\2pdata\BpodData\rbpm12\ImageNLWReverseOlfGNG';
%     'C:\2pdata\BpodData\ckrspm2\ImageNLWTrainOlfGNG';
%     'C:\2pdata\BpodData\m295\ImageNLWTrainStimGoNoGo\Jul13_2017';
%     'C:\2pdata\BpodData\m118\ImageNLWTrainStimGoNoGo\Jul13_2017';
%     'C:\2pdata\BpodData\m120\ImageNLWTrainStimGoNoGo\Jul13_2017';

%      'C:\2pdata\BpodData\rbpm12\ImageNLWTrainOlfGNG';
%      'C:\2pdata\BpodData\rbpm12\ImageNLWReverseOlfGNG';
%        'C:\2pdata\BpodData\rbpm13\ImageNLWTrainOlfGNG';
%      'C:\2pdata\BpodData\rbpm13\ImageNLWReverseOlfGNG';
%      'C:\2pdata\BpodData\ckrspm2\ImageNLWTrainOlfGNG';
%      'C:\2pdata\BpodData\ckrspm2\ImageNLWReverseOlfGNG';
%      'C:\2pdata\BpodData\rbpm14\ImageNLWTrainOlfGNG';
%        'C:\2pdata\BpodData\rbpm12\ImageNLWTrainOlfGNG';
%        'C:\2pdata\BpodData\rbpm13\ImageNLWTrainOlfGNG';
%        'C:\2pdata\BpodData\rbpm14\ImageNLWTrainOlfGNG';
%        'C:\2pdata\BpodData\m407\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m293\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m118\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m120\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m119\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m407\ImageNLWTrainStimGoNoGo';
         'C:\2pdata\BpodData\m295\ImageNLWTrainStimGoNoGo';}
       
for i = 1:numel(dirnames)
    dirname = dirnames{i}
    %%% Recurse through directory
    tic
    Convert2pRecursive( dirname );
    toc


    dos(['explorer.exe ', dirname])
    
end
% 
% dos(['C:\Users\dlab\Desktop\Fiji.app\ImageJ-win64.exe -macro ', fullfile(pwd, 'imagejmacros\openStack.ijm'), '  ', fullfile(pathname, [fname(1:end-4), '.tif'])])

%%
disp('NOW STARTING NLW_2p_cont')
tstart = tic;
addpath(genpath('C:\Users\dlab\src\OEGAnalyze'));
NLW_2p_preprocess;
toc(tstart)