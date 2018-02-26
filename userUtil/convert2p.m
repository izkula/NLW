direc = fullfile(''); %%% Set this if you will be using frequently

[fname, pathname] = uigetfile(['C:\2pdata\', direc, '*.sbx'])

disp(['Converting ', fname, ' to tif.'])
sbx2tif_optotune(fullfile(pathname, fname(1:end-4)));

dos(['explorer.exe ', pathname])

dos(['C:\Users\dlab\Desktop\Fiji.app\ImageJ-win64.exe -macro ', fullfile(pwd, 'imagejmacros\openStack.ijm'), '  ', fullfile(pathname, [fname(1:end-4), '.tif'])])
