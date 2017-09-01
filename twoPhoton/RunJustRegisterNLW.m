function RunJustRegisterNLW(vidpath, outdir, templatepath, imagejPath,...
                imagejMacroPath, doDebug, isContinuousTrial, startFrame, numFrames, ...
                downsample_val)
% Call ImageJ for denoising and image registration. 
 
%%% Directory to save intermediate videos (denoised, and then registered)
% outdir = fullfile(viddir, 'reg')
if ~exist(outdir, 'dir'); mkdir(outdir); end
if ~exist('doDebug', 'var') || isempty(doDebug); doDebug = false; end
if ~exist('isContinuousTrial', 'var') || isempty(isContinuousTrial); isContinuousTrial = false; end
if ~exist('startFrame', 'var') || isempty(startFrame); startFrame = 1; end
if ~exist('numFrames', 'var') || isempty(numFrames); numFrames = []; end


% % % % %%% Testing code for turboreg, originally developed by Sam
% % % % %    cmd = [imagejPath, ' -macro ./turboreg.ijm ''/home/izkula/src/OEGAnalyze/matlab/WholeCortex/twoPhoton/denoised_small.tif''']
% % % % 
% % % % %%% Not using headless mode, seems to work the most consistently
% % % % %  cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, '''']
% % % % %%%%cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
% % % % 
% % % % %%% To use xvfb instead of headless mode. May yield large speedup. 
% % % % % cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, '''']
         

useMOCO = true;
useHeadless = true
disp(templatepath)


if doDebug
    useHeadless = false
end

if useMOCO
    if isContinuousTrial
        if isempty(templatepath)
            templatepath = '0.0' %%%% Will generate the template.
            if ~useHeadless
%                 cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW_cont.ijm ''',vidpath,' ', outdir,  ' ', templatepath, ' ', num2str(startFrame), ' ', num2str(numFrames) '''']
                   cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW_cont.ijm ',vidpath,',', outdir,  ',', templatepath, ',', num2str(startFrame), ',', num2str(numFrames), ',', num2str(downsample_val) '']
            else
%                 cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW_cont.ijm ''',vidpath,' ', outdir,  ' ', templatepath, ' ', num2str(startFrame), ' ', num2str(numFrames) '''']
                  cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW_cont.ijm ',vidpath,',', outdir,  ',', templatepath, ',', num2str(startFrame), ',', num2str(numFrames), ',', num2str(downsample_val) '']

            end
        else
            if ~useHeadless
%                 cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate_cont.ijm ''',vidpath,' ', outdir, ' ', templatepath, ' ', num2str(startFrame), ' ', num2str(numFrames) '''']
                cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate_cont.ijm ',vidpath,',', outdir, ',', templatepath, ',', num2str(startFrame), ',', num2str(numFrames), ',', num2str(downsample_val) '']
            else
%                 cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate_cont.ijm ''',vidpath,' ', outdir, ' ', templatepath, ' ', num2str(startFrame), ' ', num2str(numFrames) '''']
                cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate_cont.ijm ',vidpath,',', outdir, ',', templatepath, ',', num2str(startFrame), ',', num2str(numFrames), ',', num2str(downsample_val) '']
            end
        end
    else
        if isempty(templatepath)
            if ~useHeadless
                cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
            else
                cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLW.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
            end
        else
            if ~useHeadless
                cmd = [imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
            else
                cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/moco_NLWwithTemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
            end
        end
    end
else
    if isempty(templatepath)
        if ~useHeadless
            cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg_noDenoiseNLW.ijm ''',vidpath,' ', outdir, ' ', '''']
        else
            cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg_noDenoiseNLW.ijm ''',vidpath,' ', outdir, ' ', '''']; %%% Uncomment this instead for headless mode
        end
    else
        if ~useHeadless
            cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg_noDenoiseNLWwithTemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
        else
            cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg_noDenoiseNLWwithTemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']

    %         cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg_noDenoiseNLWwithTemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']; %%% Uncomment this instead for headless mode
        end
    end
end
%%%%t = getCurrentTask(); t.ID
disp(outdir);


%%% Align all trials to one template image
% cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg_sessiontemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath, '''']

%   cmd = [imagejPath, ' --console -macro ', imagejMacroPath, '/turboreg.ijm ', '''/home/izkula/src/OEGAnalyze/matlab/WholeCortex/imagejMacros/denoised_small.tif''']
   
%%% For testing denoise
%cmd = [imagejPath, ' --headless --console -macro ./process2p.ijm ''folder=../folder1 parameters=a.properties output=../samples/Output'' ']
%   cmd = [imagejPath, ' --console -macro ', imagejMacroPath, '/process2p.ijm ''folder=../folder1 parameters=a.properties output=../samples/Output'' ']

disp(cmd)
system(cmd);
disp(['done: ', outdir]);
