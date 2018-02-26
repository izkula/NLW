function [] = grin2p_AppendSessions(basePath, out_path)

f = dir(fullfile(basePath,'*','*'));

%find unique mice
subjects = get_unique_subjects(f);

%open the images for each 


for ii = 1:numel(subjects) %for each mouse
    subjectName = subjects{ii};
    
for i = 1:numel(f) %for each session in the 2p data folder
    sessionDir = fullfile(f(i).folder, f(i).name)
  
        if regexp(subjectName, sessionDir(1:4));
            image_path = fullfile(sessionDir,'z0_processed')
            
            cmd = ['~/Fiji.app/ImageJ-linux64', ' -macro', ' ~/src/NLW/imagejMacros/grin_2p_append.ijm ', ...
            image_path ',' out_path];
        
            
    
    

        end   
    end
end

    
    
    



end
