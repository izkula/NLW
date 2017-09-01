fnames=dir;
fnames=fnames(3:end-2);
filenames={};
for a=1:length(fnames)
    filenames{a}=fnames(a).name;
end

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/tiffiles/20150523/';

Miji;
%filename_tif=input('filename: ');

%for a=1:length(filenames)
for a=1:length(filenames)
    filename_tif=filenames{a};
    
    tiff_info=imfinfo(filename_tif);
    nframes=length(tiff_info);
    W=tiff_info(1).Width;
    H=tiff_info(1).Height;
    
    %h1=waitbar(0,'Opening tif');
    mov=zeros(H,W,nframes);
    parfor b=1:nframes
        mov(:,:,b)=imread(filename_tif,b);
        %waitbar(b/nframes,h1);
    end
    movm=squeeze(mean(mov,3));
    
    savename=filename_tif;
    savename=strcat(pathname,savename);
    savename=strrep(savename,'.tif','_mc.tif');
    
    MIJ.createImage('target',mov,true);
    MIJ.createImage('mean',movm,true);
    MIJ.run('Single Registration v2 8 3',['choose=[' savename ,'] choose=[' pathname ']']);
    MIJ.run('Close All');
end