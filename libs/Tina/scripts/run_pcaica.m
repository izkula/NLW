function run_pcaica(pathname,filedate,fn,flims,nPCs,dsamp,badframes,mu,doSubtractMean)

if ~exist('doSubtractMean', 'var') || isempty(doSubtractMean)
    doSubtractMean = false;
end

currpath=strcat(pathname,filedate);

% create output directory if it doesn't exist
outputdir=strcat(currpath,'/timecourses/');
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

% Open tif file and perform first batch PCA
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(strcat(currpath,'/tiffiles/',fn,'.tif'), flims, nPCs, dsamp, outputdir, badframes);
[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(strcat(currpath,fn,'.tif'), flims, nPCs, dsamp, outputdir, badframes, doSubtractMean);


% Choose which PC's to use
[PCuse] = CellsortChoosePCs(strcat(currpath,fn,'.tif'), mixedfilters);

% Set parameters for ICA
nIC=length(PCuse);
ica_A_guess=randn(length(PCuse), nIC);
termtol=1e-6;
maxrounds=1000;

% Perform ICA
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);

save(strcat(outputdir,fn,'_ICA.mat'),'PCuse','ica_sig','ica_filters','ica_A','numiter','movm')
