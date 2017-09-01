function h = OverlayBoundaries(im, BW)


 [B, L, N] = bwboundaries(BW);
 h = imagesc(im); hold on;
 for k=1:length(B),
    boundary = B{k};
    if(k > N)
        plot(boundary(:,2),...
            boundary(:,1),'g','LineWidth',1);
    else
        plot(boundary(:,2),...
            boundary(:,1),'r','LineWidth',1);
    end
end