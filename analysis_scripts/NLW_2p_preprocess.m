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
doExtractTifs = 0

if doExtractTifs
    %convert from sbx to tiffs
    convert2pFolder()
    
    %remove SBX files
    RemoveSBXFiles(fullfile(basePath))   
end

%% Crop, DownSample and Register Images

doGRINImages = 0

if doGRINImages
    %downsample data spatially, crop, and register within a session
    grin_2p_CropReduceRegisterMoco(basePath)
end

%% Register and Aggregate all the sessions together
doAggregateSessions = 0 %manual right now due to dumbness

if doAggregateSessions
    %downsample data spatially, crop, and register within a session
    grin_2p_AggregateSessions(basePath)
end

%%

f = dir(fullfile(basePath,'*','*'));

for i = 1:numel(f)
    sessionDir = fullfile(f(i).folder, f(i).name)
    image_path = [sessionDir '/z0_processed.tif']
    out_path = [sessionDir '/z0_proc_reduced.tif']
    
    %check that tif existd
    if exist([f(i).folder '/' f(i).name '/' 'z0_processed.tif'])  && ~exist(out_path)
    %Run macro to load, downsample and register image stack
    cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_reduce.ijm ', ...
        image_path ',' out_path]

    fprintf(image_path)
    system(cmd)
     end
        
    %run CNMFE    
    if ~exist([f(i).folder '/' f(i).name '/z0_proc_reduced_source_extraction'],'dir')
        try
            demo_large_data_2p([f(i).folder '/' f(i).name '/' 'z0_proc_reduced.tif']);
        end
    end    
        
end
%% Upsample Traces



%% Append Traces Together



%% Organize Metadata




