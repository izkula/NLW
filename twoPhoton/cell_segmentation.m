%%%%% USER INPUTS %%%
clear all;

% alpha = .1; beta = .7; gamma = 0.2; %USER

doLoadFrames = false;
path_name = '~/Dropbox/forSethu/vglut_1215_5';
% path_name = '~/Dropbox/forSethu/gad2_1220_7';

plane_id = 1;

tic
if doLoadFrames
    %%% Load Data
    N_skip = 1;
    n_planes = 1;
   
    image_name = 'registered49.tif';
    imagefile_directory = [path_name '/' image_name];
    im_orig = imread(imagefile_directory,1); % USER file location of the image
    N_frames = 1;
    while(imread(imagefile_directory,N_frames))
        N_frames = N_frames+1;
    end
    N_frames = N_frames - 1;
    f = im_orig; %%% reads first frame here, size of f is 512 by 512
    imsize_x = size(f,1);
    imsize_y = size(f,2);
    f_base = 5000 * ones(imsize_x,imsize_y);
    begin_plane = 1;
    clear f;
    plane_id = 1;
    im_orig = double(imread(imagefile_directory,plane_id));
    %im_orig = im2double(im_orig);

    %%% Compute video statistics
    f(:,:,plane_id) = im_orig;
    f_min(1:imsize_x,1:imsize_y,plane_id) = im_orig;
    f_sum(1:imsize_x,1:imsize_y,plane_id) = im_orig;
    f_sum2(1:imsize_x,1:imsize_y,plane_id) = im_orig.^2;
    f_sum3(1:imsize_x,1:imsize_y,plane_id) = im_orig.^3;
    slice_cnt = 1;
    for i= 2:N_frames   
        ftemp = imread(imagefile_directory,i); % USER
        ftemp = double(ftemp);
        %f_pixel1(slice_cnt) = ftemp(276,376);
        %f_pixel2(slice_cnt) = ftemp(293,264);
        slice_cnt = slice_cnt + 1;
        f_sum(:,:,plane_id) = f_sum(:,:,plane_id) + ftemp;
        f_sum2(:,:,plane_id) = f_sum2(:,:,plane_id) + ftemp.^2;
        f(:,:,plane_id) = max(f(:,:,plane_id), ftemp);  %%% new f contains the max pixel intensity of frames 1 to n
        f_min(:,:,plane_id) = min(f_min(:,:,plane_id), ftemp);
        f_skew(:,:,plane_id) = f_sum2(:,:,plane_id) + ftemp.^3;
    end
    f_sum = f_sum / slice_cnt;
    f_sum2 = f_sum2 / slice_cnt;
    f_skew = f_skew / slice_cnt;


    %%% preliminary filter, not necessary
    plane_id
    f_mean = mean(mean(f(:,:,plane_id))); %%% spatial average of f, f_mean is a scalar
    f_max = max(max(f(:,:,plane_id))); %%% this is the maximum intensity


    for i = 1:imsize_x
        for j = 1:imsize_y        
            f_mu = f_sum(i,j,plane_id);
            f_std = sqrt(double( f_sum2(i,j,plane_id) - f_mu^2));
            f_stdev(i,j) = f_std;
            f_skewness = f_skew(i,j,plane_id) - 3 * f_mu * f_std^2 - f_mu^3;
            f_skew(i,j) = f_skewness;
        end
    end
  
else

    f_skew = double(imread(fullfile(path_name, 'maxSkew00001.tif')));
    f_stdev = sqrt(double(imread(fullfile(path_name, 'maxVar00001.tif'))));
    f_mu = double(imread(fullfile(path_name, 'maxMean00001.tif')));
    f_max = double(imread(fullfile(path_name, 'maxProj00001.tif')));   
    imsize_x = size(f_skew,1);
    imsize_y = size(f_skew,2);
end
figure; imshow(f_skew,[]); % displays skewness image
figure; imshow(f_stdev,[]); % displays stdev image
toc

%% Merge stats images
tic
param_threshold = 0.85; %USER
std_cutoff = mean(mean(f_stdev))+2*std(f_stdev(:)); % USER
skew_cutoff = mean(mean(f_skew))-2*std(f_skew(:)); % USER
max_cutoff = mean(mean(f_max)) - 2*std(f_max(:))*param_threshold;

radius = 4;
max_threshed = imopen(adaptivethresh(f_max, 10, -20, 'gaussian', 'relative'), strel('disk', radius));
stdev_threshed = imopen(adaptivethresh(f_stdev, 10, -20, 'gaussian', 'relative'), strel('disk', radius));
skew_threshed = imopen(adaptivethresh(f_skew, 10, -20, 'gaussian', 'relative'), strel('disk', radius));

% f_new = alpha*max_threshed ...
%         + beta*stdev_threshed ...
%         + gamma*skew_threshed;

% f_new = stdev_threshed.*skew_threshed.*max_threshed;
% f_new = (max_threshed + stdev_threshed.*skew_threshed) >= 1;
f_new = (max_threshed.*stdev_threshed + stdev_threshed.*skew_threshed) >= 1;

figure; imshow(f_new,[]); % displays new filtered image
figure, imagesc(f_new.*f_max +0.8*f_max)

%% Perform segmentation
[bwd_cells,N_distinct_cells] = bwlabel(f_new);
% cc = bwconncomp(f_new);
centroids = zeros(N_distinct_cells,2);
for i = 1:N_distinct_cells
   ss =  regionprops(bwd_cells == i);
   centroids(i,:) = ss.Centroid;
end
figure;imagesc(bwd_cells); % displays reconstructed segmentation
figure; imagesc(f_max); % displays max intensity
hold on; plot(centroids(:,1), centroids(:,2), 'ro')



orig_bwd_cells = bwd_cells;
%% Adjust segmentation
close all
bwd_cells = orig_bwd_cells;
cellPositions = select_cells(f_max, f_skew, num2cell(centroids,2));
[new_points, removed_points, new_inds, removed_inds] = GetPointChanges(cellPositions', num2cell(centroids,2));
%%% remove those points from the labeled image
for i = 1:numel(removed_inds)
    bwd_cells(bwd_cells == removed_inds(i)) = 0; 
end
%%% add new points to the labeled image
new_radius = 10;
for i = 1:numel(new_inds)
    new_cells = zeros(size(bwd_cells));
    new_cells(round(new_points(i,2)), round(new_points(i,1))) = 1;
    new_cells = imdilate(new_cells, strel('disk', new_radius));
    bwd_cells = bwd_cells + new_inds(i)*new_cells;
end
figure, imagesc(bwd_cells)


%% Perform segmentation
%{
%%%% simple edge detection %%%

f_local = f_new; %% f now contains the new image
toc
h = d2gauss(imsize_x, imsize_y, 2.0); % 3.0 is the standard deviation

%imshow(filter2(h,f),[]);
new_image_stage1 = filter2(h,f_local);

new_image_stage2 = filter2(h,new_image_stage1); %% to be explored, how many filters is optimal?

new_image_stage3 = edge(new_image_stage2,'canny',0.2); % to be explored, is 0.2 the best? 0.15/0.25?

%%% now perform cell segmentation

h = d2gauss(imsize_x,imsize_y,1.0); % 3.0 is the standard deviation, this filter is only used for the purpose of visualization
BWdfill = imfill(filter2(h,new_image_stage3), 'holes');
%%%%%%%%%%%%%%% use bwdlabel to find distinct cells inside an image
%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
im_orig = imread(imagefile_directory,plane_id); % USER file location of the image
f_post(:,:,1) = im_orig; %%% reads first frame here, size of f is 512 by 512
imsize_x = size(f_post,1);
imsize_y = size(f_post,2);
f_base = 5000 * ones(imsize_x,imsize_y);
cnt = 2;

for i=1:length(BWdfill)
    for j=1:length(BWdfill)
        if (BWdfill(i,j) > 0.2)
            new_var(i,j,plane_id)=1;
        else
            new_var(i,j,plane_id)=0;
        end
    end
end
[bwd_cells{plane_id},N_distinct_cells(plane_id)] = bwlabel(new_var(:,:,plane_id));
figure;imagesc(bwd_cells{plane_id}); % displays reconstructed segmentation
figure; imshow(f,[]); % displays max intensity
%}