function run_findcellmasks(pathname,filedate,fn,smwidth,thresh,arealims,plotting)

currpath=strcat(pathname,filedate);

load(strcat(currpath,'/timecourses/',fn,'_ICA.mat'),'movm','ica_filters');
outputdir=strcat(currpath,'/timecourses/');

masktype=0;
% Use these lines for semi-automatic ROI selection
[ica_segments, segmentlabel, segcentroid, ica_centers,icfilter_cell] = CellsortSegmentation_semiauto(ica_filters, smwidth, thresh, arealims, plotting, movm, masktype);

if masktype==0
    cellmask=ica_segments;
else
    cellmask=ica_centers;
end


% %Use for manual ROI drawing
% [ica_segments, segmentlabel, cellmask, segcentroid] = CellsortSegmentation_manual(ica_filters, smwidth, thresh, arealims, plotting);
% for a=1:size(cellmask,3)
%     newcellmask(a,:,:)=cellmask(:,:,a);
% end
% cellmask=newcellmask;

% Plot cells one by one
tiff_info=imfinfo(strcat(currpath,fn,'.tif'));
nframes=length(tiff_info);
W=tiff_info(1).Width;
H=tiff_info(1).Height;
npix=W*H;

save(strcat(outputdir,fn,'_cellmasks_og.mat'),'ica_segments','ica_centers','segmentlabel','cellmask','segcentroid','movm','icfilter_cell')

% % % figure
% % % colormap(gray)
% % % imagesc(squeeze(movm));
% % % hold on
% % % colorcode=rand(size(cellmask,1),3);
% % % for a=1:size(cellmask,1)
% % %     contour(squeeze(cellmask(a,:,:)),1,'color',colorcode(a,:));
% % %     axis([0 W 0 H])
% % %     hold on
% % %     text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','m','color','r')
% % % end

