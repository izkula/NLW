function [] = grin2pReduceImageSize(basePath)
%
f = dir(fullfile(basePath,'*','*'));


for i = 1:numel(f)
    sessionDir = fullfile(f(i).folder, f(i).name)
    
    
    %check that tif existd
    if exist([f(i).folder '/' f(i).name '/' 'z0'],'dir') && ~exist([f(i).folder '/' f(i).name '/' 'z0_processed.tif'])
        
    image_path = [sessionDir '/z0']
    out_path = [sessionDir '/z0_processed.tif']
    
     
    %Run macro to load, downsample and register image stack
    cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_process.ijm ', ...
        image_path ',' out_path]

    fprintf(image_path)
    system(cmd)
    fprintf(out_path)
    end
end

