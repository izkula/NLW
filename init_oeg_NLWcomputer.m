%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath L R A magentaC greenC blueC orangeC
% basePath = 'C:\2pdata\'
basePath = 'C:\'
% bpodImagePath = '/media/Data/data/BpodImageData'
bpodImagePath = 'C:\2pdata\BpodData'
% bpodDataPath = '~/Dropbox/Bpod_r0_5/Data';
bpodDataPath = 'C:\Users\dlab\Dropbox\Bpod\Data';
imagejPath = 'C:\Users\dlab\Desktop\Fiji.app\ImageJ-win64.exe'
%imagejPath = '/home/izkula/Software/Fiji.app/ImageJ-linux64'
resultsPath = '~\Dropbox\NLW_results'
imagejMacroPath = 'C:\Users\dlab\src\OEGAnalyze\matlab\WholeCortex\imagejMacros'

% temp = load('registration/selectionPoints.mat'); 
temp = load('registration\selectionPoints_barrel_LR.mat'); 
L = temp.L;
try
    R = temp.R;
catch
    R = [];
end
A = load('registration\atlas.mat');

magentaC = [1 0 1];
greenC = [0 1 0];
blueC = [0 0 1];
orangeC = [1 .6 0];