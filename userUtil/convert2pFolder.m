direc = fullfile(''); %%% Set this if you will be using frequently

dirname = uigetdir(['C:\2pdata\', direc])


%%% Recurse through directory
tic
Convert2pRecursive( dirname );
toc


%dos(['explorer.exe ', dirname])
