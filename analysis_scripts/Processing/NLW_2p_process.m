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

% Crop, DownSample and Register Images

doGRINImages = 1

if doGRINImages
    %downsample data spatially, crop, and register within a session
    grin_2p_CropReduceRegisterMoco(basePath)
end

%% Then run CNMFE

%% Organize meta data
grin_2p_OrganizeMetaData;

%% extract movement data
%to complete, because some of the 8arms don't have this make this of
%secondary importance
%grin_2p_ExtractMovement;


%% put it all together in a data struct

