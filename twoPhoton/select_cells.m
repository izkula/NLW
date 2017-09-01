function valid_points = select_cells(im_maxproj, movie, existing_points, labeledROIs)
% Usage: 
%   left click: add a new point, or drag to reposition an existing point
%   right click: delete a point
%   any keyboard key press: play movie
%   left click out-of-points: exit
% 
% Args:
%   im_maxproj: maximum project image, dimension R x S array
%   movie: dimension R x S x T array
%   existing_points: (optional) cell array containing 1x2 double arrays of existing
%     point coordinates
% 
% Returns:
%   valid_points: cell array containing 1x2 double arrays of selected
%     coordinates
    
if nargin < 3
    existing_points = {};
end
if nargin < 2
    movie = []
end
if nargin < 1
    % if no inputs given, make up some demo data
    N = 512;
    im_maxproj = 2*rand(N,N);
    movie = rand(N,N,20);
end
%%
f = figure('Position', [1000, 1000, 1300, 1300]);
% h = imshow(im_maxproj,[]);
% h = imagesc(im_maxproj); colormap('gray')
h = OverlayBoundaries(im_maxproj, labeledROIs); 
set(h,'cdata', im_maxproj);
title('when done, click outside the image');

% pre-populate existing points
points = {}; num_points = 0;
for i = 1:length(existing_points)
    pe = existing_points{i};
    if (pe(1) > 0)
        p = impoint(gca, pe(1), pe(2));
%     else
%         p = impoint(gca, 0, 0);
%         delete(p);
    else
        p = {};
    end
    num_points = num_points + 1;
    points{num_points} = p;            
end

% continue until user specifies otherwise
for i = 1:500000
    w = waitforbuttonpress;
    if w == 0 % mouse button clicked        
        cursor_pos=get (gca, 'CurrentPoint');
        % if cursor is left or above image, then end point selection
        if cursor_pos(1,1) <= 0 || cursor_pos(1,2) <= 0 || cursor_pos(1,1) > size(im_maxproj,1) || cursor_pos(1,2) > size(im_maxproj,2)
            % only return points that haven't been deleted
            valid_points = {};  k = 0;
            for j = 1:length(points)            
% % %                 if isvalid(points{j}) % check if point hasn't been deleted
% % %                     k = k+1;
% % %                     valid_points{k} = points{j}.getPosition;                    
% % %                 end
                if ~isempty(points{j}) && isvalid(points{j}) % check if point hasn't been deleted
                    valid_points{j} = points{j}.getPosition;   
                    k = k+1;
                else
                    valid_points{j} = [-1,-1];
                end
            end
            disp(['Selected ' num2str(k) ' points.']);
            close(f)
            return % valid_points
        end
                
        % check if a point exists at that cursor location
        point_exists = false;
        for j = 1:length(points)            
            if ~isempty(points{j}) && isvalid(points{j})
                pe = points{j}.getPosition;
                t = 3; % minimum allowed space between points
                if abs(cursor_pos(1,1) - pe(1)) < t && abs(cursor_pos(1,2) - pe(2)) < t                    
                    point_exists = true;
                    % if right click and point exists, delete it
                    if strcmp(get(f,'SelectionType'),'alt')
                        disp('Removing point!')
                          delete(points{j});
                    end
                end
            end
        end
        % add a point if not already existing and left click
        if ~point_exists && strcmp(get(f,'SelectionType'),'normal')
            p = impoint(gca, cursor_pos(1,1), cursor_pos(1,2));            
            num_points = num_points + 1;
            p.setColor('green');
            points{num_points} = p;            
        end        
    else % a keyboard key pressed
        % play movie
        duration = 0.7; % seconds
        for k = 1:size(movie,3)
            set(h,'cdata', squeeze(movie(:,:,k))/max(max(movie(:,:,1)))*max(max(im_maxproj(:))));
%             set(h, 'cdatamapping', 'scaled')
%             set(h, 'caxis', [min(min(movie(:,:,1))), max(max(movie(:,:,1)))]);

            pause(duration / size(movie,3));
        end
        % show max project image again
        set(h,'cdata', im_maxproj);
%         set(h, 'caxis', [min(im_maxproj(:)), max(im_maxproj(:))]);
    end
end