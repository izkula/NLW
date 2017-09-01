function status = ExtractTifs(tifPath, targetPath, firstNumber, doCrop)
% Extract a multipage tif into individual tifs and append into a target
% folder. Start the numbering at firstNumber.
% Requires imagemagick to be installed and accessible by command line. 
% http://www.imagemagick.org/script/index.php

if ~exist('doCrop', 'var') || isempty(doCrop)
    doCrop = false
end

%  convert input -crop WxH+X+Y +repage output
 
cmd = ['convert -scene ', num2str(firstNumber), ' ', tifPath, ' ', targetPath, '/%07d.tif']

status = system(cmd)

if status == 127
    error('Could not find `convert` command for use in ExtractTifs. Make sure that imagemagick is installed')
end
%%%% Be sure to hang until the process is completed. 

% convert -scene 10 '/home/izkula/Data/data/BpodData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session2/Trial00001/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session2_Trial00001/z0/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session2_Trial00001__001_z0.tif' /home/izkula/Data/data/BpodData/m594/ImageNLWTrainStimGoNoGo/May04_2017/Session2/Trial00001/m594_ImageNLWTrainStimGoNoGo_May04_2017_Session2_Trial00001/z0/single%d.tif
