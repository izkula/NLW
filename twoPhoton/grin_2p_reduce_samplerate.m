function [] = grin_2p_reduce_samplerate(resultsPath,subject)
%UNTITLED Summary of this function goes here

if ~exist(fullfile(resultsPath, subject, 'z_cat_red.tif'),'file')
    if exist(fullfile(resultsPath, subject, 'z_cat_full.tif'),'file')
        
       image_path = fullfile(resultsPath, subject, 'z_cat_full.tif')
       out_path = fullfile(resultsPath, subject, 'z_cat_red.tif')
       
       cmd = ['~/Fiji.app/ImageJ-linux64','-macro','~/src/NLW/imagejMacros/grin_2p_reduce_sampelrate.ijm',...
           image_path ',' out_path]
       
       fprintf(image_path)
       system(cmd)
       fprintf(out_path)
        
        
    else
        
       fprintf('The full image stack has not been created')
    
    
    
    end
    



end





end

