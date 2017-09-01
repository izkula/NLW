function processedDataDir = GetProcessedDataDir2p(origDataPath, basePath, bpodImagePath)
%%% Generates pathname for saving out processed data (so that it is 
%%% held in a separate location than the original data, but with an
%%% equivalent pathname).

% global basePath bpodImagePath

processedDataDir = strrep(origDataPath, bpodImagePath,  ...
                          fullfile(basePath, 'processedData', 'BpodImageData'));