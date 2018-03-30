%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath trial_data frameTicks timeTicks frameSampleRate
basePath = '~/2pdata/'

%bpodImagePath = 'C:\2pdata\BpodData'
%bpodDataPath = 'C:\Users\dlab\Dropbox\Bpod\Data';
imagejPath = '~\Fiji.app\ImageJ-linux64.exe'
resultsPath = '~/2pdata/Results'
imagejMacroPath = fullfile(pwd, 'imagejMacros')


% LOG OF GRIN DATA FOR USEABLE DATASETS %
%object = 1; social = 2; toymouse = 3; blank = 4;
trial_data = {...
    'm873' '02/19/2018' [1 2 1 4 2 1 2 3]
    'm872' '02/19/2018' [1 2 1 4 2 1 2 3]
    'm876' '02/18/2018' [1 2 1 4 2 1 2 3]
 };

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