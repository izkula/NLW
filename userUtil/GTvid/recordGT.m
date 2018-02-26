function [ data, time ] = recordGT( vid, nframes, fps, savename)
%RECORDGT Record video (keep in memory)
%   nframes - number of frames
%   fps - frames per second
%   savename - if [] then saves to memory and returns in 'data', else saves
%              to disk

%%% Note: use imaqtool if you want a prebuilt gui!

if ~exist('savename', 'var') || isempty(savename)
    saveToDisk = false;
else
    saveToDisk = true;
end

stop(vid);
src = getselectedsource(vid);
vid.ReturnedColorspace = 'grayscale';
vid.FramesPerTrigger = nframes;
src.AcquisitionFrameRateAbs = fps;

if saveToDisk
    vid.LoggingMode = 'disk';
     logfile = VideoWriter([savename, '.avi'], 'Grayscale AVI');
%     logfile = VideoWriter([savename, '.mj2'], 'Archival');

    vid.DiskLogger = logfile;
else
    vid.LoggingMode = 'memory';
end

if ~saveToDisk
    preview(vid)
end

tic
start(vid);


if saveToDisk
    wait(vid, nframes/fps+5);
    while(vid.FramesAcquired ~= vid.DiskLoggerFrameCount)
        pause(0.1);
    end
    disp('Done saving')
    data = [];
    time = [];
    toc
else
    while(vid.FramesAvailable < nframes)
         vid.FramesAvailable
    end
    toc
    [data, time] = getdata(vid, nframes);
    stoppreview(vid);
end




