%%LOADS THE NIDAQ OUTPUT AND EXTRACTS THE RUNNING SIGNAL, DOWNSAMPLES TO 31
%%HZ AND MOVES IT TO THE RESULTS DIRECTORY

init_samz

f = dir(fullfile(basePath,'*','*','*'))

filepaths = {}; j = 1
for i = 1:numel(f)
    if contains(f(i).name,'nidaq.mat')
        filepaths{j} = fullfile(f(i).folder, f(i).name);
        j = j+1;

    end
end

clear f

%for each nidaq file
for i = 1:numel(filepaths)
    f = filepaths{i}
    try
    
    [t, ch] = LoadNidaqOutput(f,1);
    
    nCh = numel(ch)
    
    catch
        fprintf('No Nidaq file')
    end
    
    
end



    