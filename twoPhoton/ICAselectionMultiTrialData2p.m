function ICAselectionMultiTrialData2p(dataName, imagePath, varargin)
% Computes ICA components used for later cell segmentation.
%
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath

%%% Function options
p = inputParser();
p.addParameter('zplanes', [], @isnumeric);
p.addParameter('flims', []);
p.addParameter('nPCs', [], @isnumeric);
p.addParameter('dsamp', [], @isnumeric);
p.addParameter('badframes', [], @isnumeric);
p.addParameter('mu', [], @isnumeric);

p.parse(varargin{:});

zplanes = p.Results.zplanes;
flims = p.Results.flims;
nPCs = p.Results.nPCs;
dsamp = p.Results.dsamp;
badframes = p.Results.badframes;
mu = p.Results.mu;


processedDataDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath)

if ~exist(processedDataDir, 'dir'); mkdir(processedDataDir); end


nzplanes = 1;
if ~isempty(zplanes)
    nzplanes = numel(zplanes)
end

for zz = 1:nzplanes

    if ~isempty(zplanes)
        zplane = zplanes(zz);
        pathname = fullfile(processedDataDir, dataName, ['z', num2str(zplane)],'/'); 
    else
        zplane = [];
        pathname = fullfile(processedDataDir, dataName, ['z0/']); 
    end
    
    filefolder = 'reg/';
    filename = 'registered_noDenoise';
    
    doSubtractMean = true;
%     pathname = '/home/izkula/Data/processedData/BpodImageData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session1/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session1/z0/'
%     filedate = 'reg_meansub/'
%     filename = 'registered_noDenoise_meansub'
    
     run_pcaica(pathname,filefolder,filename,flims, nPCs,dsamp,badframes,mu, doSubtractMean);
end

end

