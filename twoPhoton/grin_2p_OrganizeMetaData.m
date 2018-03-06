function [] = grin_2p_OrganizeMetaData(basePath)
%
f = dir(fullfile(basePath,'*','*'));





slice_counter = 0
frames = [];
datasets = {};

for i = 1:numel(f)
    if strmatch(f(i).name(1),'m')

        sessionDir = fullfile(f(i).folder, f(i).name)
    
    %load metadata
    list_dir = dir(sessionDir)
    
    for i = 1:numel(list_dir)
        if ~isempty(strmatch(list_dir(i).name(1), 'm')) && list_dir(i).bytes < 2000 %if it starts with m and is a small file
            meta_data_fname = list_dir(i).name
            load(fullfile(sessionDir, meta_data_fname))
            break
        end
    end
      
    
    %get the number of slices in the processed, registed tiff stack
    I = imfinfo(fullfile(sessionDir,'z0_proc_reg.tif'));
    nSlices = numel(I);
    
    %add this to slice counter
    slice_counter = slice_counter + nSlice
    
    %calculate the 
    
    
    
    
    
    
    
    end
end
 