function traces = ExtractTraces2p(cellMasks, vid)
% cellMasks is a labeled image of cell ROIs
% Extract a matrix of traces from the provided video. 
%
% Returns: matrix where each row corresponds to a cell, each column is a
%          time point

Ncells = max(cellMasks(:));
traces = zeros(Ncells, size(vid,3));

for kk =1:Ncells
%     kk
    inds = find(cellMasks == kk);
    Ninds = numel(inds);
    
    trace = zeros(size(vid, 3),1);
    for i = 1:Ninds
        [x1,x2] = ind2sub(size(cellMasks), inds(i));
        trace = trace + squeeze(vid(x1, x2, :)); %%% I manually checked this ordering. 
    end
    traces(kk, :) = trace/Ninds;
end

function SubtractNeuropil(images)
disp('subtract neuropil');

function neuropilMasks = FindNeuropilMasks(image)
disp('find neuropil mask');



% function traces = ExtractTraces2p(cellPositions, vid, cellRadius)
% % Extract a matrix of traces from the provided video. 
% %
% % Returns: matrix where each row corresponds to a cell, each column is a
% %          time point
% 
% cp = round(cellPositions); 
% 
% h = fspecial('disk', cellRadius);
% v = imfilter( vid, h, 'same'); 
% 
% traces = zeros(size(cp, 1), size(vid,3));
% for kk =1:size(cp, 1)
%     traces(kk,:) = v(cp(kk, 2), cp(kk, 1), :); % I double checked this ordering. 
% end