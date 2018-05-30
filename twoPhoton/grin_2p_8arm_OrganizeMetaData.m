init_samz
%
f = dir(fullfile(basePath,'*','*'));

%%
names = {}; j = 1

for i = 1:numel(f)
    try
        %get concatenated file names
        if contains(f(i).name,'concatenated')
            names{j} = f(i).name(1:4);
            j = j+1;
        end
    end
end
%%

for iii = 1:numel(names) %for all mice
    
    %gather all the filenames  for this mouse
    filepaths = {}; j = 1;
    for i = 1:numel(f)
        if ~contains(f(i).name, 'contatenated')
        if strmatch(names{iii}, f(i).name)
            filepaths{j} = fullfile(f(i).folder, f(i).name);
           j = j+1 ; %counter
        end
    end
    end

    %Define and reset important variables of metadata
    slice_counter = 0 %number of slices in the whole image stack
    frames = []; %frame numbers for trial starts and ends
    datasets = {}; 
    for ii = 1:numel(filepaths)

        %load metadata
        list_dir = dir(filepaths{ii});

        for i = 1:numel(list_dir)
            if ~isempty(strmatch(list_dir(i).name(1), 'm')) && list_dir(i).bytes < 2000 %if it starts with m and is a small file
                meta_data_fname = list_dir(i).name
                load(fullfile(filepaths{ii}, meta_data_fname))
                break
            end
        end


        %get the number of slices in the processed, registed tiff stack
        I = imfinfo(fullfile(filepaths{ii},'z0_processed.tif'));
        nSlices = numel(I);

        frames = [frames; info.frame + slice_counter];
        
        plot(frames)
        
        %add this to slice counter
        slice_counter = slice_counter + nSlices
        
        data{ii} = filepaths{ii};

    end
        

meta.frames = frames;

if ~exist(fullfile('~', '2presults', [names{iii} '8arm']))
    mkdir(fullfile('~', '2presults', [names{iii} '_8arm']))
end

save(fullfile('~', '2presults',[names{iii} '_8arm'],'metadata.mat'),'meta')

end


    
