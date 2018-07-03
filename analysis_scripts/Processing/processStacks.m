%%% Process stacks captured by Neurolabware microscope

addpath(genpath('~/src/NLW'))


%%% Define stack locations
dataDir = '/media/Data/data/'
saveDir = '/media/Data/processedData/'

do_all_in_dir = true

if ~do_all_in_dir
    %%% {'date', 'stack_name', frames_per_z, overall_depth};
    data = {...
    %     {'20180529', 'm194_x1500_y1500_z-583_stack920_1-4x', 30, 300};
        {'20180529', 'm194_x1500_y1636_z-158_stack920_3-4x_120um', 30, 120};
        {'20180529', 'm194_x1500_y-1500_z-418_stack920_1-4x', 30, 300};
    }
else
    date = '20180529'
    frames_per_z = 30
    dirlist = dir(fullfile(dataDir, date));
    data = {}
    for d = 1:numel(dirlist)
        if strfind(dirlist(d).name, 'stack')
            if strfind(dirlist(d).name, '120um')
                depth = 120
            else
                depth = 300
            end
            data{end+1} = {date, dirlist(d).name, frames_per_z, depth};
        end
    end
end

%% Load stacks

do_register = false

for i = 1:numel(data)
    d = data{i} 
    frames_per_z = d{3}

    imageDir = fullfile(dataDir, d{1}, d{2}, 'z0')
    currImages = LoadImageStack(imageDir);

    if do_register
       disp('Not yet implemented')
    end

    nz = floor(size(currImages, 3)/frames_per_z);
    stack = zeros(size(currImages, 1), size(currImages, 2), nz);
    for z = 1:nz
       stack(:,:,z) = mean(currImages(:,:,(z-1)*frames_per_z+1:z*frames_per_z), 3);
    end
    stack = stack(100:500, 100:500,:);

    savePath = fullfile(saveDir, d{1}, d{2})
    if ~exist(savePath, 'dir')
       mkdir(savePath)
    end
    saveFile = fullfile(savePath, 'stack.mat')
    save(saveFile, 'stack')
   
%    MakeAVI(stack, 10, fullfile(savePath, 'stack.avi', cmap, maxval)
    close all
    MakeAVI(stack, 10, fullfile(savePath, 'stack.avi'))

    SaveTiffStack(stack, fullfile(savePath, 'stack'))
end