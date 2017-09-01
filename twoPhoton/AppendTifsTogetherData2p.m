function AppendTifsTogetherData2p(dataName, imagePath, varargin)
% Motion corrects and registers 2p data across trials within a bpod 
% behavioral session.
%
% Saves registered videos to processedData folder. 
%
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session2' 
% imagePath: i.e.
%   basePath/BpodImageData/vglut1m15/Image2PShapeOlfGNG/Dec10_2015/Session2
%%%

global basePath bpodImagePath imagejPath imagejMacroPath

%%% Function options
p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('zplanes', [], @isnumeric);
p.parse(varargin{:});

numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
zplanes = p.Results.zplanes;

processedDataDir = GetProcessedDataDir2p(imagePath, basePath, bpodImagePath)

if ~exist(processedDataDir, 'dir'); mkdir(processedDataDir); end

%%% Load data
if ~isempty(numTrials)
    subSeq = 1:numTrials
else
    subSeq = []
end


nzplanes = 1;
if ~isempty(zplanes)
    nzplanes = numel(zplanes)
end

for zz = 1:nzplanes

    if ~isempty(zplanes)
        zplane = zplanes(zz);
    else
        zplane = [];
    end
    
    targetDir = fullfile(imagePath, dataName, ['z', num2str(zplane)]);
    
    %%% If metadata file exists here, then this was recorded continuously
    %%% and you do not actually need to do any appending. You just need to
    %%% copy the trialStart times from the metadata file. 
    metaFileName = fullfile(imagePath, dataName, [dataName, '__001.mat']);
    if exist(metaFileName, 'file')
       actuallyDoAppend = false;
       mfile = load(metaFileName);
       trialStartInds = mfile.info.frame(mfile.info.event_id == 1);
       trialEndInds = [trialStartInds(2:end)-1; trialStartInds(end) + 1];
       %%% Double check that didn't seem to skip any trials.
       trialTimes = diff(trialStartInds);
       trialTimes = trialTimes(1:end-1);
       if max(trialTimes) - median(trialTimes) > median(trialTimes)
           warning(['DROPPED a trial start ind in ', dataName]);
           didNLWRecordStartFrame = zeros(size(trialStartInds));
       else
           didNLWRecordStartFrame = ones(size(trialStartInds));
       end
       
    else
       actuallyDoAppend = true;
    end
    
    
    if actuallyDoAppend
        if exist(targetDir, 'dir')
           rmdir(targetDir, 's');
        end
        mkdir(targetDir)
    
        [~, vidpaths, viddirs] = Load2pImageTrials(imagePath, dataName, ...
                                                   'subSeq', subSeq, ...
                                                   'doLoad', false, ...
                                                   'isNLW', isNLW, ...
                                                   'zplane', zplane); % Obtain the paths to the video stacks

        disp(['Num vid dirs: ', num2str(numel(viddirs))])

        trialStartInds = zeros(numel(vidpaths), 1);
        trialEndInds = zeros(numel(vidpaths), 1);
        didNLWRecordStartFrame = zeros(numel(vidpaths), 1);


        for kk = 1:numel(vidpaths)
            disp(['Trial ', num2str(kk)])
            if isNLW
                fsplit = strsplit(vidpaths{kk}, '/');
                metaDataFileName = fullfile('/', fsplit{1:end-2}, [fsplit{end-2}, '__001.mat']);
                m = load(metaDataFileName);
                if isfield(m.info, 'frame')
                    trialStart = m.info.frame(1);
                    didNLWRecordStartFrame(kk) = 1;
                else
                    trialStart = 0;
                end
            else
               trialStart = 0; 
            end

            firstNumber = GetFinalImageNumberInDir(targetDir) + 1
            tifPath = vidpaths{kk}
            ExtractTifs(tifPath, targetDir, firstNumber)
            %%%% Need to pause until the process is completed?
            lastNumber = GetFinalImageNumberInDir(targetDir)

            trialStartInds(kk) = firstNumber + trialStart;
            trialEndInds(kk) = lastNumber;
        end
    end
    
    save([targetDir, '_trialStarts.mat'], 'trialStartInds', 'trialEndInds', 'didNLWRecordStartFrame')

end

end

