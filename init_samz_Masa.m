%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath trial_data frameTicks timeTicks frameSampleRate
basePath = '~/2pdata/Masa/'

%bpodImagePath = 'C:\2pdata\BpodData'
%bpodDataPath = 'C:\Users\dlab\Dropbox\Bpod\Data';
imagejPath = '~\Fiji.app\ImageJ-linux64.exe'
resultsPath = '~/Dropbox/2p_Masa/Results/'
imagejMacroPath = fullfile(pwd, 'imagejMacros')


% LOG OF GRIN DATA FOR USEABLE DATASETS %
%object = 1; social = 2; toymouse = 3; blank = 4;
trial_data = {...
    '6517Het' 'Spontaneous' 'GCaMP6f'
 };

%for 8 arm
frameSampleRate = 31;

timeTicks = frameTicks / frameSampleRate;