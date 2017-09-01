function [new_points, removed_points, new_inds, removed_inds] = GetPointChanges(newPoints, origPoints)
% Given two cell arrays, each containing entries of x, y point coordinates
% returns a vectors of points (each row is a point)
% that have been added to origPoints and
% points that have been removed from origPoints to yield newPoints. 
% removed_inds indicates the index into origPoints of the points that have
% been removed. 

NP = cell2mat(newPoints);
OP = cell2mat(origPoints);

if isempty(OP)
    new_points = NP;
    new_inds = 1:size(NP,1);
    removed_points = [];
    removed_inds = [];
else
%     [new_points, new_inds] = setdiff(NP, OP, 'rows');
%     [removed_points, removed_inds] = setdiff(OP, NP, 'rows');
    new_points = NP(size(OP, 1)+1: size(NP,1),:)
    new_inds = size(OP,1)+1:size(NP,1)
    
    removed_points = [];
    removed_inds = [];
    for row = 1:size(OP, 1)
        if NP(row, 1) < 0 && OP(row, 1) > 0
            removed_points = [removed_points; OP(row,:)]
            removed_inds = [removed_inds; row]
        end
    end
end



