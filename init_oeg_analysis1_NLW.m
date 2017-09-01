%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath L R A magentaC greenC blueC orangeC
basePath = '/media/Data/'
% bpodImagePath = '/media/Data/data/BpodImageData'
bpodImagePath = '/media/Data/data/BpodData'
% bpodDataPath = '~/Dropbox/Bpod_r0_5/Data';
bpodDataPath = '~/Dropbox/Bpod/Data';
imagejPath = '/home/izkula/Software/ImageJ/ImageJ'
%imagejPath = '/home/izkula/Software/Fiji.app/ImageJ-linux64'
resultsPath = '~/Dropbox/NLW_results'
imagejMacroPath = '/home/izkula/src/OEGAnalyze/matlab/WholeCortex/imagejMacros'

% temp = load('registration/selectionPoints.mat'); 
temp = load('registration/selectionPoints_barrel_LR.mat'); 
L = temp.L;
try
    R = temp.R;
catch
    R = [];
end
A = load('registration/atlas.mat');

magentaC = [1 0 1];
greenC = [0 1 0];
blueC = [0 0 1];
orangeC = [1 .6 0];