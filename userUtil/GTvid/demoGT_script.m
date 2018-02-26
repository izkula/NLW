%%%% Capture video with the GT2750 gige camera


%% Initialize
vid = initGT();

%% Preview

preview(vid);

%% Capture video

fps = 10;
nframes = 1000;

savePath = [];
savePath = 'cux2_fork_2'

[data, time] = recordGT(vid, nframes, fps, savePath);

data = squeeze(data(:,:,1,:));


%% Clean up
delete(vid); clear vid




