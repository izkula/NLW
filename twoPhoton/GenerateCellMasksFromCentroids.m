function labeledMasks = GenerateCellMasksFromCentroids(centroids, cellRadius, labels, baseImage)
%%% Add labeled cell masks corresponding to 'centroids' location to an
%%% existing image of cell masks. 

%%% centroids - a vector where each row is the coordinates of a point
%%% cellRadius - (int) default radius of a cell outline
%%% labels - (vector of ints) number corresponding to each of the new centroids
%%% baseImage - an image: either zeros, or contains already labeled cell masks to
%%% which the new masks will be added (labels should not overlap with the
%%% existing mask labels).

labeledMasks = baseImage;
if numel(labels > 0)
    for i = 1:numel(labels)
        masks = zeros(size(baseImage));
        if centroids(i,1) > 0 && centroids(i, 2) > 0
            masks(round(centroids(i,2)), round(centroids(i,1))) = 1;
            masks = imdilate(masks, strel('disk', cellRadius));
            labeledMasks = labeledMasks + labels(i)*masks;
        end
    end
end