function [labeledROIs, cellCentroids] = AutoSegmentCells(f_max, f_var, f_skew, thresh)
%%% Segments cells (from two-photon videos) based on three image
%%% statistics:
%%%     - the maximum intensity projection of each pixel across all time points
%%%     - the variance of each pixel, maximum projected across all
%%%     behavioral trials in a session
%%%     - the skewness of each pixel, maximum projected across all
%%%     behavioral trials in a session
%%%
%%% Outputs:
%%%     - labeledROIs - an image, the same size as the input images with
%%%     masks each labeled corresponding to each cell
%%%     - the centroid positions of each labeled cell, a vector where each
%%%     row is the coordinates (correct order for accessing into an image
%%%     is  image(round(centroids(i,2)), round(centroids(i,1))) 

f_stdev = sqrt(f_var);

%%% Merge stats images
if ~exist('thresh', 'var') || isempty(thresh)
    thresh = -20; %USER
end
% std_cutoff = mean(mean(f_stdev))+2*std(f_stdev(:)); % USER
% skew_cutoff = mean(mean(f_skew))-2*std(f_skew(:)); % USER
% max_cutoff = mean(mean(f_max)) - 2*std(f_max(:))*thresh;

radius = 4;
max_threshed = imopen(adaptivethresh(f_max, 10, thresh, 'gaussian', 'relative'), strel('disk', radius));
stdev_threshed = imopen(adaptivethresh(f_stdev, 10, thresh, 'gaussian', 'relative'), strel('disk', radius));
skew_threshed = imopen(adaptivethresh(f_skew, 10, thresh, 'gaussian', 'relative'), strel('disk', radius));


f_new = (max_threshed.*stdev_threshed + stdev_threshed.*skew_threshed) >= 1;

figure; imshow(f_new,[]); % displays new filtered image
figure, imagesc(f_new.*f_max +0.8*f_max)


%%% Perform segmentation
[bwd_cells,N_distinct_cells] = bwlabel(f_new);
centroids = zeros(N_distinct_cells,2);
for i = 1:N_distinct_cells
   ss =  regionprops(bwd_cells == i);
   centroids(i,:) = ss.Centroid;
end

cellCentroids = centroids;
labeledROIs = bwd_cells;

% figure;imagesc(bwd_cells); % displays reconstructed segmentation
% figure; imagesc(f_max); % displays max intensity
% hold on; plot(centroids(:,1), centroids(:,2), 'ro')


