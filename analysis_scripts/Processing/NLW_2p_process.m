%%%% Preprocess data on samz computer before transferring to data
%%%% computer. Only modify this file on the NLW computer. 

 clear all
 cd('~/src/NLW/')
addpath(genpath('.')); 
colordef white; 
init_samz;

%%% globals are generated using computer dependent init file (i.e. initOEG or init_oeg_analysis1_NLW)
global basePath  imagejMacroPath imagejPath resultsPath

%% Extract TIFFs
doExtractTifs = 1

if doExtractTifs
    %convert from sbx to tiffs
    convert2pFolder()
    
    %remove SBX files
    RemoveSBXFiles(fullfile(basePath))   
end
%
% Crop, DownSample and Register Images

doGRINImages = 1

if doGRINImages
    %downsample data spatially, crop, and register within a session
    grin_2p_CropReduceRegisterMoco(basePath,'20180522')
end

%% MAKE AVERAGE PROCESSED VIDEOS


%% Then run CNMFE

%% Organize meta data and extract movement data

processMetaData;


grin_2p_8arm_OrganizeMetaData;

%grin_2p_running_OrganizeMetaData;

%then move contours, neuron, processed and traces to results folder



%% put it all together in a data struct

