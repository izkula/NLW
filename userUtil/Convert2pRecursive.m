function Convert2pRecursive( dirname )
%Recursively convert all .sbx files to tifs in a
%directory tree.

D = dir(dirname);

for i = 3:numel(D) %%par
    d = D(i);
    if d.isdir
        Convert2pRecursive(fullfile(d.folder, d.name));
    else
        if ~isempty(strfind(d.name, '.sbx'))
            tifExists = 0;
            for kk = 1:numel(D) %%% Check that it hasn't already been converted
                if strcmp('z0', D(kk).name)
                    tifExists = 1;
                end
            end
            if~tifExists
                disp(['Converting ', fullfile(d.folder, d.name(1:end-4))]);
                try
                    sbx2tif_optotune(d.folder, d.name(1:end-4));
                catch
                    disp(['Could not convert ', d.folder, d.name]);
                end
            end
        end
    end
end


end

