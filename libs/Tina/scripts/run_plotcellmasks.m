function run_plotcellmasks(pathname,filedate,fn)

currpath=strcat(pathname,filedate);
load(strcat(currpath,'/timecourses/',fn,'_cellmasks.mat'),'movm','cellmask','segcentroid');

% Plot cells
tiff_info=imfinfo(strcat(currpath,fn,'.tif'));
nframes=length(tiff_info);
W=tiff_info(1).Width;
H=tiff_info(1).Height;
npix=W*H;


labels = zeros(W, H);
for a = 1:size(cellmask, 1)
    labels(cellmask(a, :,:) > 0) = a+100;
end
figure, imagesc(labels)

% figure;
% colormap(gray)
% imagesc(squeeze(movm));
% hold on
% colorcode=rand(size(cellmask,1),3);
% for a=1:size(cellmask,1)
% %     contour(squeeze(cellmask(a,:,:)),1,'color',colorcode(a,:),'linewidth',2);
%     axis([0 W 0 H])
%     hold on
%     contour(squeeze(cellmask(a,:,:)),[1 1],'color',colorcode(a,:));
%     text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% end
% 
% title(strcat(num2str(size(cellmask,1)), 'Cells'));


