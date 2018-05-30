function [] = grin_2p_CropReduceRegisterMoco(basePath)
%
f = dir(fullfile(basePath,'*','*'));

z = {'z0' 'z1' 'z2'};
for i = 1:numel(f)
    sessionDir = fullfile(f(i).folder, f(i).name)
    
    for j = 1:3
    
    %%%%%%%check that z0 tif existd%%%%%%
    if exist([f(i).folder '/' f(i).name '/' z{j}],'dir') && ~exist([f(i).folder '/' f(i).name '/' z{j} '_processed.tif'])
        
    f_images = dir([f(i).folder '/' f(i).name '/' z{j}])  
    nImages = size(f_images, 1)-2
        
    image_path = [sessionDir '/' z{j}]
    out_path = [sessionDir '/' z{j} '_processed.tif']
    
     
    %Run macro to load, downsample and register image stack
    cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_process.ijm ', ...
        image_path ',' out_path ',' nImages]

    fprintf(image_path)
    system(cmd)
    fprintf(out_path)
    
    end
    
 
    
    end
end

