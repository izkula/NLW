function [RGB_xy] = overlayTargets(im, cellpos, filter_image, regions, label)
% Overlay numbered labels on image of cells, to indicate which trace corresponds
% to which cell. 
% cellpos is a vector where each row is cell coordinates. 

x = cellpos(:, 1);
y = cellpos(:, 2);
gamma =  1; %1/1.8;
if nargin < 2
    filter_image = false;
end
if nargin < 3
    color_cells = {'Yellow'};
else
   color_cells =  {};
   for i = 1:size(regions,1)
       if mod(regions(i),5) == 1
           color_cells{i} = 'cyan';
       elseif mod(regions(i),5) == 2
           color_cells{i} = 'red';
       elseif mod(regions(i),5) == 3
           color_cells{i} = 'blue';
       elseif mod(regions(i),5) == 4
           color_cells{i} = 'green';
       elseif mod(regions(i),5) == 0
           color_cells{i} = 'yellow';
       end
   end
end
if nargin < 5
    label = [1:size(x,1)];
end

%% median filter
% dim = 1;
% if filter_image
%     w = 7; % size of median filter kernel
%     for i = 1:size(im,3)
%         im(:,:,i) = medfilt2(im(:,:,i),[w w]);
%     end
%     maxproj_xy = max(im, [], 3);
%     minproj_xy = min(im, [], 3);
%     maxproj_xy = maxproj_xy - minproj_xy;
%     %maxproj_xy = max(0,maxproj_xy - prctile(maxproj_xy(end/4:end*3/4,end/4:end*3/4),5));
%     
%     maxproj_yz = squeeze(max(im, [], dim))';
%     minproj_yz = squeeze(min(im, [], dim))';
%     maxproj_yz = maxproj_yz - minproj_yz;
%     
%     figure, imshow(maxproj_xy,[]);
%     rg = [prctile(maxproj_xy(:),15) prctile(maxproj_xy(:),99.99)];
%     ima = imadjust(maxproj_xy/max(maxproj_xy(:)),rg,[0 1], 1/2.2);
%     figure, imshow(ima);
% 
% else
    maxproj_xy = im;
    rg = [0 1];

% end
% %% fix aspect ratio to be isotropic
% maxproj_yz = imresize(maxproj_yz, [round(size(maxproj_yz,2)*200/(2048*6.5/16*35/50) ) size(maxproj_yz,2)]);
% z = round(z* size(maxproj_yz,1)/51);

%%
r = 5;
position = [x, y, repmat(r, size(x))];

scale = 2;
RGB_xy = insertObjectAnnotation(uint8(255*imadjust(imresize((double(maxproj_xy)/max(maxproj_xy(:))), scale*size(maxproj_xy), 'nearest'),rg, [0 1], gamma)), ...
    'circle', position*scale, label, 'TextColor', 'black','Color', color_cells);

% figure(1), close(1), figure(1); 
% imshow(RGB_xy, []); colormap jet;