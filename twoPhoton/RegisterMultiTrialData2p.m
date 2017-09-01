function RegisterMultiTrialData2p(dataName, imagePath, varargin)
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
p.addParameter('doMotionCorrect', true, @islogical);
p.addParameter('doTrialRegister', true, @islogical);
p.addParameter('doDenoise', false, @islogical);
p.addParameter('doRegister', true, @islogical);
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('useTemplate', false, @islogical); %%% Register to template (median of first trial in session)
p.addParameter('zplanes', [], @isnumeric);
p.addParameter('doDebug', false, @islogical); 
p.addParameter('isContinuousTrial',false,@islogical);
p.addParameter('startFrame',1,@isnumeric);
p.addParameter('numFrames',[],@isnumeric);
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});

doMotionCorrect = p.Results.doMotionCorrect;
doTrialRegister = p.Results.doTrialRegister;
doDenoise = p.Results.doDenoise;
doRegister = p.Results.doRegister;
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
useTemplate = p.Results.useTemplate;
zplanes = p.Results.zplanes;
doDebug = p.Results.doDebug;
isContinuousTrial = p.Results.isContinuousTrial;
startFrame = p.Results.startFrame;
numFrames = p.Results.numFrames;
tDownsampleFactor = p.Results.tDownsampleFactor;

% processedDataDir = strrep(imagePath, bpodImagePath,  ...
%                           fullfile(basePath, 'processedData', 'BpodImageData'))
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
    
    if ~isContinuousTrial
        [~, vidpaths, viddirs] = Load2pImageTrials(imagePath, dataName, ...
                                                   'subSeq', subSeq, ...
                                                   'doLoad', false, ...
                                                   'isNLW', isNLW, ...
                                                   'zplane', zplane); % Obtain the paths to the video stacks
    else
        viddirs = {fullfile(imagePath, dataName, ['z', num2str(zplane)])};
        vidpaths = {fullfile(viddirs{1}, '0000001.tif')};
    end

    disp(['Num vid dirs: ', num2str(numel(viddirs))])

    % template = median(vids{1}, 3);template_dir = fullfile(viddirs{1}, 'template');
    % mkdir(template_dir)
    % template_path =  fullfile(template_dir, 'template');
    % template_fname = SaveTiffStack( template, template_path); %%% This didnt actually seem to work once loaded into imageJ. Probably a file format issue. Better to just generate median image in imagej even if it is slightly slower. 

    if ~isNLW
        for jj = 1:numel(viddirs)
            FilterImageDirectory(viddirs{jj}, 'Ch2');
        end
    end

    template_path = vidpaths{1};
    sessStartTime = tic;

    downsample_val = (1/tDownsampleFactor)*nzplanes

    %%% Generate template image
    if useTemplate
        currdir = viddirs{1};
        
        if isempty(zplanes)
            outdir = GetProcessedDataDir2p(fullfile(imagePath, 'template'), basePath, bpodImagePath);
        else
            outdir = GetProcessedDataDir2p(fullfile(imagePath, 'template', ['z', num2str(zplane)]), basePath, bpodImagePath);
        end
        mkdir(outdir);
        if ~exist(fullfile(outdir, 'median.tif'), 'file')
            templateStartFrame = 1;
            templateNumFrames = 200;
            RunJustRegisterNLW(template_path, outdir, [], imagejPath, imagejMacroPath, ...
                doDebug, isContinuousTrial, templateStartFrame, templateNumFrames, downsample_val);
            template_path = fullfile(outdir, 'median.tif');
        else
            template_path = fullfile(outdir, 'median.tif');
        end
    end



    %%% Copy over small metadata mat file to processed data dir (so everything
    %%% is in the same place).
    % if isNLW
    %     for jj = 1:numel(viddirs)
    %         currdir = viddirs{jj}
    %         metaoutdir = GetProcessedDataDir2p(currdir, basePath, bpodImagePath);
    %         if ~exist(metaoutdir, 'dir')
    %             mkdir(metaoutdir);
    %         end
    %         metadatafile = dir(fullfile(currdir, '*.mat'));
    %         copyfile(fullfile(currdir, metadatafile.name), fullfile(metaoutdir, metadatafile.name))
    %     end
    % end


    if doDenoise %%% Not used with NLW currently
        if ~doRegister
            disp('running denoise')
            tic
            for jj = 1:numel(viddirs)
                currdir = viddirs{jj};
                outdir = GetProcessedDataDir2p(fullfile(currdir, 'reg'), basePath, bpodImagePath);
                if strcmp(vidpaths{jj}(end-3:end), '.tif')
                    RunDenoise(vidpaths{jj}, outdir, template_path, imagejPath, imagejMacroPath);
                else
                    disp(['Did NOT denoise: ', vidpaths{jj}]);
                end
            end
            toc
        end

        disp('running registration')
        tic
        for jj = 1:numel(viddirs)
            currdir = viddirs{jj};
            outdir = GetProcessedDataDir2p(fullfile(currdir, 'reg'), basePath, bpodImagePath);            
            RunRegister(vidpaths{jj}, outdir, template_path, imagejPath, imagejMacroPath);
            disp(num2str(toc(sessStartTime)));
        end
        toc

    else
        disp('registering without denoise')
    %     parfor jj = 1:numel(viddirs)
        for jj = 1:numel(viddirs)
            currdir = viddirs{jj}
            datestr(now)
            outdir = GetProcessedDataDir2p(fullfile(currdir, 'reg'), basePath, bpodImagePath)
            tic
            
            if strcmp(vidpaths{jj}(end-3:end), '.tif')
                if ~exist(fullfile(outdir, 'registered_noDenoise.tif'),'file')
                    if isNLW
                        if ~useTemplate
                            currtemplate_path = [];
                        else
                            currtemplate_path = template_path;
                        end
                        
                        if isContinuousTrial
                            if isempty(numFrames)
                               numFrames = GetFinalImageNumberInDir(viddirs{jj});
                            end
                        end
                        
                        RunJustRegisterNLW(vidpaths{jj}, outdir,  currtemplate_path, ...
                                                 imagejPath, imagejMacroPath, doDebug, ...
                                                 isContinuousTrial, startFrame, numFrames, ...
                                                 downsample_val);  
                    else
                        RunJustRegister(vidpaths{jj}, outdir,  template_path, ...
                                                 imagejPath, imagejMacroPath);  
                    end
                end
                disp(num2str(toc(sessStartTime)));
                

                %%% Copy some files over to the processed data dir. 
                if ~isempty(zplanes)
                    currdir = currdir(1:end-3); %%% Remove the zlabel (assumes fewer than 10 planes...should be reasonable....)
                end
                metaoutdir = GetProcessedDataDir2p(currdir, basePath, bpodImagePath);
                metadatafile = dir(fullfile(currdir, '*1.mat'));
                if ~isempty(metadatafile)
                    try
                        copyfile(fullfile(currdir, metadatafile.name), fullfile(metaoutdir, metadatafile.name))
                    catch
                        disp(['Could not copy ', fullfile(currdir, metadatafile.name)])
                    end
                end
                
%                 metadatafile = dir(fullfile(currdir, '*trialStarts.mat'));
                metadatafile = dir(fullfile(currdir, ['z', num2str(zplane), '_trialStarts.mat']));
                if ~isempty(metadatafile)
                    try
                        copyfile(fullfile(currdir, metadatafile.name), fullfile(metaoutdir, metadatafile.name))
                    catch
                        disp(['Could not copy ', fullfile(currdir, metadatafile.name)])
                    end
                end
            else
                disp(['Did NOT register: ', vidpaths{jj}]);
            end
        end
        disp('elapsed time');
        toc
    end

end

end

