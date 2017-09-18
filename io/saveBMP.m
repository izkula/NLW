function [] = saveBMP(I, cmap, savepath, cmin, cmax)
%%% cmap should be something like hsv(256) or fire(256)
 

if ~exist('cmax', 'var') || isempty(cmax)

    cmax = max(I(:));

end

if ~exist('cmin', 'var') || isempty(cmin)

    cmin = min(I(:));

end

    

I = (I-cmin)./(cmax - cmin);

I = min(1, max(0, I));

Irecolored = zeros(size(I,1), size(I,2), 3);

Irecolored(:,:,1) = reshape(cmap(round(I(:).*255)+1, 1), size(I));

Irecolored(:,:,2) = reshape(cmap(round(I(:).*255)+1, 2), size(I));

Irecolored(:,:,3) = reshape(cmap(round(I(:).*255)+1, 3), size(I));

 

imwrite(Irecolored, savepath);