%%% Define global variables (for data from Neurolabware 2p)

global basePath bpodDataPath bpodImagePath imagejPath imagejMacroPath resultsPath
% basePath = 'C:\2pdata\'
basePath = '/home/sam/2pNLW';
bpodImagePath = [basePath '/BpodData/Data'];
bpodDataPath = '~/Dropbox/Bpod/Data';
imagejPath = '/opt/Fiji.app/ImageJ-linux64.exe';
resultsPath = '~/Dropbox/NLW_results';
% imagejMacroPath ='C:\Users\dlab\src\OEGAnalyze\matlab\WholeCortex\imagejMacros'
imagejMacroPath = fullfile(pwd, 'imagejMacros');
