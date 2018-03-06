function [] = grin_2p_AggregateSessions(basePath)
%UNTITLED Summary of this function goes here %%THIS DOES NOT WORK%%

f = dir(fullfile(basePath, '*', '*'));


for i = 1:numel(f)
    
   
    sessionDir = fullfile(f(i).folder, f(i).name)
    image_path = [sessionDir '/z0_processed.tif']
    out_path = [sessionDir '/z0_processed_reg.tif']

    if exist(image_path)

    if  find(regexpi(f(i).name, '001'))
        %make the registration directory
        if ~exist([basePath '/Results/' f(i).name(1:4) '/Registration/' ])
            mkdir([basePath '/Results/' f(i).name(1:4) '/Registration/'])
        end
        
        out_path_template = (['/home/svesuna/2pdata/Results/' f(i).name(1:4) '/Registration/template.tif'])
        %Run macro to load the first session and make a template image
        cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_makeTemplate.ijm ', ...
        image_path ',' out_path_template]

    system(cmd)
    fprintf(out_path)
        
    end
    
    %then register each session to this template image and save out
    
    %Run macro to load, downsample and register image stack
    cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_registerToTemplate.ijm ', ...
        image_path ',' out_path_template ',' out_path]

    fprintf(image_path)
    system(cmd)
    fprintf(out_path)
    end
end

end

