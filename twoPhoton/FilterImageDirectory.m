function FilterImageDirectory(viddir, template)

    trials = dir(viddir);
    
    garbageDir = fullfile(viddir, 'garbage');
    mkdir(garbageDir);
    
    for i = 1:numel(trials)
        name = trials(i).name;
        if ~isdir(fullfile(viddir, name))
            if isempty(strfind(name, template)) && isempty(strfind(name, '.xml'))  && isempty(strfind(name, '.env'))
                movefile(fullfile(viddir, name), fullfile(garbageDir, name));
            end
        end
    end

end

