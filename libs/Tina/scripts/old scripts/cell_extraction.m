%% PCA/ICA

% Set various parameters for opening tif file
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/20150518/';
fn_base='wtnac_reward002';

outputdir=strcat(pathname,'timecourses/');
flims=[];
nPCs=20;
dsamp=[];
badframes=[];

% Open tif file and perform first batch PCA
[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_RS(strcat(pathname,'tiffiles/',fn), flims, nPCs, dsamp, outputdir, badframes);

% Choose which PC's to use
[PCuse] = CellsortChoosePCs_RS(fn, mixedfilters);

% Set parameters for ICA
mu=0.1; % keep between 0.1-0.2
nIC=length(PCuse);
ica_A_guess=randn(length(PCuse), nIC);
termtol=1e-6;
maxrounds=1000;

% Perform ICA
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA_RS(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);

save(strcat(outputdir,strrep(fn,'.tif','_ICA')),'PCuse','ica_sig','ica_filters','ica_A','numiter')

%% Cell Segmentation

close all

% Set parameters for smoothing
smwidth=0;
thresh=10; % set high to avoid overlapping cells
arealims=50;
plotting=1;

% Use these lines for automatic ROI selection
[ica_segments, segmentlabel, segcentroid] = CellsortSegmentation_RS(ica_filters, smwidth, thresh, arealims, plotting);
cellmask=ica_segments;

% % Make ROI's smaller
% for a=1:size(cellmask,1);
%     currcellmask=cellmask(a,:,:);
%     currcellmask(currcellmask<11)=0;
%     cellmask(a,:,:)=currcellmask;
% end

% Plot cells one by one
tiff_info=imfinfo(strcat(pathname,'tiffiles/',fn));
nframes=length(tiff_info);
W=tiff_info(1).Width;
H=tiff_info(1).Height;
npix=W*H;

figure
colormap(gray)
imagesc(squeeze(movm));
hold on
colorcode=rand(size(cellmask,1),3);
for a=1:size(cellmask,1)
    contour(squeeze(cellmask(a,:,:)),1,'color',colorcode(a,:));
    axis([0 W 0 H])
    hold on
    text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','m','color','r')
end

del=[];
del=input('any cells to delete? ');
cellmask(del,:,:)=[];
segcentroid(del,:)=[];

save(strcat(outputdir,strrep(fn,'.tif','_cellmasks')),'ica_segments','segmentlabel','cellmask','segcentroid','movm')