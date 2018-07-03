%%%% Preprocess data on samz computer before transferring to data
%%%% computer. Only modify this file on the NLW computer. 

 clear all
 cd('~/src/NLW/')
addpath(genpath('.')); 
colordef white; 
init_EmilyS;

%%% globals are generated using computer dependent init file (i.e. initOEG or init_oeg_analysis1_NLW)
global basePath  imagejMacroPath imagejPath resultsPath

%% Extract TIFFs
doExtractTifs = 0
if doExtractTifs
    %convert from sbx to tiffs
    convert2pFolder()
    
    %remove SBX files
    RemoveSBXFiles(fullfile(basePath))   
end
%
% Crop, DownSample and Register Images
doprocessGRINStackinImageJ = 1


if doprocessGRINStackinImageJ
    %downsample data spatially, crop, and register within a session
    grin_2p_CropReduceRegisterMoco(basePath,'BpodData')
end

%%RUN CNMFE 

%% Organize meta data and extract movement data

%processMetaData_habenula

%grin_2p_running_OrganizeMetaData;

%then move contours, neuron, processed and traces to results folder



%% put it all together in a data struct

