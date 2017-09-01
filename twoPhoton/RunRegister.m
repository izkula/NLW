function RunRegister(vidpath, outdir, templatepath, imagejPath, imagejMacroPath)
% Call ImageJ for denoising and image registration. 
 
%%% Directory to save intermediate videos (denoised, and then registered)
% outdir = fullfile(viddir, 'reg')
if ~exist(outdir, 'dir'); mkdir(outdir); end

%%% Testing code for turboreg, originally developed by Sam
%    cmd = [imagejPath, ' -macro ./turboreg.ijm ''/home/izkula/src/OEGAnalyze/matlab/WholeCortex/twoPhoton/denoised_small.tif''']

%%% Not using headless mode, seems to work the most consistently
%  cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, '''']
%%%%cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']

%%% To use xvfb instead of headless mode. May yield large speedup. 
% cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, '''']

doSilent = false;
if doSilent
    cmd = ['xvfb-run -a ', imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, ' ', templatepath '''']
else
   cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg.ijm ''',vidpath,' ', outdir, ' ', templatepath ''''] 
end
disp(outdir);


%%% Align all trials to one template image
% cmd = [imagejPath, ' -macro ', imagejMacroPath, '/turboreg_sessiontemplate.ijm ''',vidpath,' ', outdir, ' ', templatepath, '''']

%   cmd = [imagejPath, ' --console -macro ', imagejMacroPath, '/turboreg.ijm ', '''/home/izkula/src/OEGAnalyze/matlab/WholeCortex/imagejMacros/denoised_small.tif''']
   
%%% For testing denoise
%cmd = [imagejPath, ' --headless --console -macro ./process2p.ijm ''folder=../folder1 parameters=a.properties output=../samples/Output'' ']
%   cmd = [imagejPath, ' --console -macro ', imagejMacroPath, '/process2p.ijm ''folder=../folder1 parameters=a.properties output=../samples/Output'' ']

system(cmd);
disp(['done: ', outdir]);
