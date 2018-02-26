function [] = RemoveSBXFiles(basePath)
%Removes the SBX file after processing to TIFFS

%
sbxfiles = dir(fullfile(basePath,'**','*.sbx'));

for i = 1:numel(sbxfiles)
    %check that tiffs exist
    if exist(sbxfiles(i).folder,'dir') && ~isempty(fullfile(sbxfiles(i).folder,'z0'))
        
    %delete sbx file
    delete(fullfile(sbxfiles(i).folder, sbxfiles(i).name))
    fprintf(['deleting ' fullfile(sbxfiles(i).folder, sbxfiles(i).name)])
    
    else
    fprintf([sbxfiles(i).name ' has not been converted to TIFFS'])
    convert2pFolder()
        
    end
end
end

