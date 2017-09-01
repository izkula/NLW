function [ output_args ] = PlotCellLocs(cellLocs, cellSubset,  labels, clim, refPoints,cmap,showHist)
if nargin < 6
    cmap = 'jet';
end
if nargin < 7
    showHist = true;
end
if isempty(cellSubset)
    cellSubset = 1:size(cellLocs,1);
end

figure('position',[100 100 520 500]);
maxX = max(cellLocs(:,1));
maxY = max(cellLocs(:,2));
cellLocs(:,1) = cellLocs(:,1)/maxX;
cellLocs(:,2) = cellLocs(:,2)/maxY;
currCellLocs = cellLocs(cellSubset,:);

subplot(4,4,[5,6,7, 9,10,11, 13,14,15]);

xl = [0 1.1]; yl = [0 1.1];
scatter(cellLocs(:,1), cellLocs(:,2), 20, [0.9 0.9 0.9], 'filled'); caxis([-0.5 0.5]); hold on;
scatter( currCellLocs(:,1), currCellLocs(:,2),20,labels,'filled'); hold on;
%text(refPoints(1,1), refPoints(1,2),'A','fontsize',16);
%text(refPoints(2,1), refPoints(2,2),'P','fontsize',16);   
%xlim([0 1]); ylim([0 1]);
caxis(clim); colormap(jet); 
colormap(cmap); axis off; axis tight; 
circle((refPoints(1,1)/maxX), (refPoints(1,2)/maxY)/2+0.05, 0.43); hold on;

plot([refPoints(1,1)/maxX refPoints(2,1)/maxX],[refPoints(1,2)/maxY*.9, refPoints(2,2)/maxY*1.8],'--k','linewidth',0.1);
xlim(xl); ylim(yl);

xc = xlim;
yc = ylim;

if showHist
    subplot(4,4,1:3);
    histogram(currCellLocs(:,1),25,'edgecolor',[0 0 0],'displaystyle','stairs'); axis off; axis tight;
    xlim(xl);    
    ylim([0 50]);

end

if showHist
    subplot(4,4,[8; 12; 16]);
    
    histogram(yc(2)-currCellLocs(:,2),25,'edgecolor',[0 0 0],'displaystyle','stairs');
    axis off;view(90,90); axis tight;   
    xlim(yl);
    ylim([0 50]);
end

end

