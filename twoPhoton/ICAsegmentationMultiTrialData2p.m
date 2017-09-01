function ICAsegmentationMultiTrialData2p(dataName, imagePath, varargin)
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
p.addParameter('smwidth', []);
p.addParameter('thresh', [], @isnumeric);
p.addParameter('arealims', [], @isnumeric);
p.addParameter('plotting', [], @isnumeric);

p.parse(varargin{:});

zplanes = p.Results.zplanes;
smwidth = p.Results.smwidth;
thresh = p.Results.thresh;
arealims = p.Results.arealims;
plotting = p.Results.plotting;


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
    
    run_findcellmasks(pathname,filefolder,filename,smwidth, thresh, arealims, plotting);

    run_deletecellmasks(pathname,filefolder,filename);
end

end

