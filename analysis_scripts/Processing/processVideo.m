clear; clc; close all
init_samz
%
%claustrum optowindow files
%f = dir(fullfile(basePath,'20180507','*','*'));

%claustrum GRIN files
f = dir(fullfile(basePath,'20180522','*','*'))

j = 1
%get concatenated file names from the session directory
for i = 1:numel(f)
    try
        if contains(f(i).name,'eye')
            fnames{j} = f(i).name;
            filepaths{j} = fullfile(f(i).folder, f(i).name);
            folderpaths{j} = f(i).folder;
            
            j = j+1;
        end
    end
end

%%

for ii = 1:numel(folderpaths)
    if ~exist(fullfile(folderpaths{ii}, 'EyeVideo'))
        try
        load(filepaths{ii})
        load(fullfile(folderpaths{ii}, 'metadata.mat'))

        data = squeeze(data(:,:,1,meta.frameStart:meta.frameEnd-1));

        options = struct('color', false,'overwrite',true,'big',true);
        tic
        saveastiff(data, [folderpaths{ii} '/EyeVideo'],options)
        toc
        catch
        end
    else
        fprintf('Already made Video')
    end
    

end
