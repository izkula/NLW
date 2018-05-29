function [cropped, x, y] = CropVid(vid, x, y)
%%% x is a list of 2 x coordinates, y is list of 2 corresponding y coords


if nargin < 2
    figure, imagesc(vid(:,:,1))
    [y,x] = ginput(2);
end

xmin = round(min(x));
xmax = round(max(x));
ymin = round(min(y));
ymax = round(max(y));

cropped = vid(xmin:xmax, ymin:ymax, :);
disp('Finished cropping')


