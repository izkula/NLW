%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath trial_data frameTicks timeTicks frameSampleRate
basePath = '~/Dropbox/2p_Claustrum_Shared/2p/'

%bpodImagePath = 'C:\2pdata\BpodData'
%bpodDataPath = 'C:\Users\dlab\Dropbox\Bpod\Data';
imagejPath = '~\Fiji.app\ImageJ-linux64.exe'
resultsPath = '~/Dropbox/2p_Claustrum_Shared/2p/Results/'
imagejMacroPath = fullfile(pwd, 'imagejMacros')


% LOG OF GRIN DATA FOR USEABLE DATASETS %
%object = 1; social = 2; toymouse = 3; blank = 4;


%for 8 arm
frameSampleRate = 31;
nFramesBaseline = 5*frameSampleRate;
nFramesAdvance = 309;
nFramesStationary = 310;
nFramesRetreat = 309;
nFramesPostTrial = 5*frameSampleRate;



%%make frame and time ticks for trial event changes
frameTicks = [nFramesBaseline nFramesBaseline+nFramesAdvance ...
    nFramesBaseline+nFramesAdvance+nFramesStationary ...
    nFramesBaseline+nFramesAdvance+nFramesStationary+nFramesRetreat]

timeTicks = frameTicks / frameSampleRate;