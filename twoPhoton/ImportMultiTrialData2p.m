function ImportMultiTrialData2p(dataName, imagePath, SessionData, varargin)
%%%
% Imports 2photon data for all trials from a given bpod behavioral session.
% Saves output traces to an hdf5 file. 
%
% Parameters:
% dataName: i.e. 'vglut1m15_Image2PShapeOlfGNG_Dec10_2015_Session1' 
% imagePath: 

%%%
   
    global basePath
    
    %%% Function options
    p = inputParser();
    p.addParameter('doMotionCorrect', false, @islogical);
    p.addParameter('doSelectCells', false, @islogical);
    
    p.parse(varargin{:});
   
    doMotionCorrect = p.Results.doMotionCorrect;
    doSelectCells = p.Results.doSelectCells;
    
    resultsDir = fullfile(basePath, 'results', dataName);
    processedDataDir = fullfile(basePath, 'processedData', dataName);    
    atlasDir = processedDataDir;    
    
    
    mkdir(resultsDir);
    mkdir(processedDataDir);


%%% Load data


%%% Register (motion correct)


%%% Select cells
% select_cells(im_maxproj, movie, existing_points)


%%% Output traces