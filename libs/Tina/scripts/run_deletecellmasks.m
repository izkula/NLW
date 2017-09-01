function run_deletecellmasks(pathname,filedate,fn)

display('Deleting cells...');

currpath=strcat(pathname,filedate);

load(strcat(currpath,'/timecourses/',fn,'_cellmasks_og.mat'));
outputdir=strcat(currpath,'/timecourses/');

% Plot cells
tiff_info=imfinfo(strcat(currpath,fn,'.tif'));
nframes=length(tiff_info);
W=tiff_info(1).Width;
H=tiff_info(1).Height;
npix=W*H;

% plot original cell masks
% % figure(100);
% % subplot(1,2,1);
% % colormap(gray)
% % imagesc(squeeze(movm));
% % hold on
% % colorcode=rand(size(cellmask,1),3);
% % for a=1:size(cellmask,1)
% %     contour(squeeze(cellmask(a,:,:)),[1 1],'color','r');
% %     axis([0 W 0 H])
% %     hold on
% %     text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','m','color','w')
% % end

% look for cellmasks that have two cells in them
for a=1:size(cellmask,1)
    currcell=logical(squeeze(cellmask(a,:,:)));
    currprops=regionprops(currcell,'area');

    if length(currprops)>1
        display(strcat('May need to split cells in ROI ',num2str(a)));
    end
end

% threshold the cell masks further
finalarea=zeros(1,size(cellmask,1));
startarea=zeros(1,size(cellmask,1));

for a=1:size(cellmask,1)
    currmask=logical(squeeze(cellmask(a,:,:)));
    temp=regionprops(currmask,'area');
    currarea=temp.Area;
    startarea(a)=currarea;
    
    thresh=0;
    while currarea>=25
        currmask=squeeze(cellmask(a,:,:));
        
        % convert cellmask values to z-scores
        temp_mean=mean(currmask(currmask>0));
        temp_std=std(currmask(currmask>0));
        currmask(currmask>0)=( currmask(currmask>0) - temp_mean ) ./ temp_std;
        
        % only keep cellmask values > thresh zscores or minimum area
        currmask(currmask<=thresh)=0;
        currmask=logical(currmask);
        
        % delete junk pixels
        currprops=regionprops(currmask,'area','pixellist');
        if length(currprops)>1
            areas = cat(1, currprops.Area);
            [~, keep]=max(areas);
            delrois=[1:length(currprops)];
            delrois(keep)=[];
            for b=1:length(delrois)
                delpix=currprops(delrois(b)).PixelList;
                for c=1:size(delpix,1)
                    currmask(delpix(c,2),delpix(c,1))=0;
                end
            end
        end
        
        % calculate new area of cell mask
        currmask=logical(currmask);
        currprops=regionprops(currmask,'area');
        currarea=currprops.Area;
        
        thresh=thresh+0.01;
    end
    
    currmask=logical(currmask);
    finalarea(a)=currarea;
    cellmask(a,:,:)=currmask;
    center=regionprops(currmask,'centroid');
    segcentroid(a,:)=center.Centroid;
end

% % % figure(200);
% % % plot(startarea,'.-');
% % % hold on
% % % plot(finalarea,'.-');
% % % legend('Starting area','Final area');
% % % xlabel('Cell number');
% % % ylabel('Area in sq. pixels');

% % % figure(300);
% % % colormap(gray)
% % % imagesc(squeeze(movm));
% % % hold on
% % % colorcode=rand(size(cellmask,1),3);
% % % for a=1:size(cellmask,1)
% % %     contour(squeeze(cellmask(a,:,:)),[1 1],'color','r');
% % %     axis([0 W 0 H])
% % %     hold on
% % %     text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','c','color','w')
% % % end


% identify cell pairs that are fewer than 5 pixels apart
distances = dist(segcentroid');
distances(linspace(1,numel(distances),length(distances))) = 100;
idx = tril(true(size(distances)), -1);
distances(idx) = 100;
[delI,delJ]=find(distances<5);

% delete the cell with the smaller area
del_all=[delI delJ];
del=[];
for a=1:size(del_all,1)
    cell1=regionprops(squeeze(cellmask(del_all(a,1),:,:)),'area');
    cell2=regionprops(squeeze(cellmask(del_all(a,2),:,:)),'area');
    areas=[cell1.Area cell2.Area];
    [tempy tempi]=min(areas);
    del(a)=del_all(a,tempi);
end
del=unique(del);
% % for a=1:length(del)
% %     contour(squeeze(cellmask(del(a),:,:)),[1 1],'color','g');
% %     axis([0 W 0 H])
% %     hold on
% % end

% del=input('Enter cellmasks to delete: ');
cellmask(del,:,:)=[];
ica_centers(del,:,:)=[];
ica_segments(del,:,:)=[];
segcentroid(del,:)=[];
icfilter_cell(del)=[];


% % figure(100);
% % subplot(1,2,2);
% % colormap(gray)
% % imagesc(squeeze(movm));
% % hold on
% % colorcode=rand(size(cellmask,1),3);
% % for a=1:size(cellmask,1)
% %     contour(squeeze(cellmask(a,:,:)),[1 1],'color','r');
% %     axis([0 W 0 H])
% %     hold on
% %     text(segcentroid(a,1), segcentroid(a,2), num2str(a), 'horizontalalignment','c', 'verticalalignment','c','color','w')
% % end


save(strcat(outputdir,fn,'_cellmasks.mat'),'ica_segments','ica_centers','segmentlabel','cellmask','segcentroid','movm','icfilter_cell')
%save(strcat(outputdir,strrep(fn,'.tif','_cellmasks')),'ica_segments','ica_centers','segmentlabel','cellmask','segcentroid','movm')
