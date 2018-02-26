function sbx2tif(folder,dname, varargin)

% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

z = sbxread(fullfile(folder, dname),1,1);
global info;

if info.volscan
%     nZplanes = numel(info.otwave_um); %%% Doesn't seem to work?
    nZplanes = 3; 
else
    nZplanes = 1;
end

if(nargin>2)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end

k = 1;
done = 0;
plane = 0; %1?

if ~exist(fullfile(folder, ['z0']), 'dir')
    while(~done && k<=N)
        try
            zfolder = fullfile(folder, ['z', num2str(plane)]);
            if ~exist(zfolder, 'dir')
                mkdir(zfolder);
            end
            
            saveAsMultiPage = false;
            if saveAsMultiPage
                outfname = fullfile(zfolder, [dname, '_z', num2str(plane)]);
            else
%               outfname = fullfile(zfolder, num2str(floor(k/nZplanes)));
                outfname = fullfile(zfolder, sprintf('%07d', (floor(k/nZplanes))));
            end

            q = sbxread(fullfile(folder, dname),k,1);
            redChan = false;
            if size(q, 1) > 1
                redChan = true;
                r = squeeze(q(2,:,:));
                zrfolder = [zfolder, '_r'];
                if ~exist(zrfolder, 'dir')
                    mkdir(zrfolder);
                    disp('saving out red channel');
                end
                outrfname = fullfile(zrfolder, sprintf('%07d', (floor(k/nZplanes))));
            end
            q = squeeze(q(1,:,:));
            
%             q = uint8(double(q)/(2^16-1))*255; %%% Convert to 8 bit (assumes filling dynamic range). 
%             q = squeeze(q(2,:,:));


            if saveAsMultiPage
                if(k<=nZplanes)
                    imwrite(q,[outfname '.tif'],'tif');
                else
                    imwrite(q,[outfname '.tif'],'tif','writemode','append');
                end
            else
                imwrite(q, [outfname, '.tif'], 'tif');
                if redChan
                    imwrite(r, [outrfname, '.tif'], 'tif'); 
                end
            end
        catch
            disp('Broke in sbx2tif_optotune')
            done = 1;
        end
        k = k+1;
        plane = mod(plane + 1, nZplanes);
    end
else
    disp(['Have already converted ', fullfile(folder, ['z0'])]);
end
    