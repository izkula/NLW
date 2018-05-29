function [pupilA, pupilR, pupilC, pupilTimes, overlayVid] = TrackPupil(vid, vidTimes, thresh, eyeMask, varargin)
%%% Tracks the location and size of the eye pupil in the provided
%%% video of the pupil, i.e. 'mantaFrames'.
%%% Returns:  normArea - area of pupil normalized to size of the eye. 
%%%           r, c - radius and center of the pupil
%%%           pupilTimes - time of each pupil frame
%%%           overlayVid - overlay of pupil circle on the video


p = inputParser();
p.addParameter('showPupilVid', false, @islogical); % Display each frame as it is processed
p.addParameter('usePupilArea', true, @islogical); % Calculate size of pupil based on it's area
p.parse(varargin{:});

showPupilVid = p.Results.showPupilVid;
usePupilArea = p.Results.usePupilArea;

warning('off', 'images:imfindcircles:warnForLargeRadiusRange');

disp(['Thresh ', num2str(thresh)]);

maxArea = sum(eyeMask(:)); % Size of the eye. 

vid = vid(:,:,1:min(numel(vidTimes), size(vid,3)));
vidTimes = vidTimes(:,1:min(numel(vidTimes), size(vid,3)));

nT = size(vid, 3);
r = zeros(nT, 1); c = zeros(nT, 2); areas = zeros(nT, 1); 

mInds = 1:nT;
overlayVid = zeros(size(vid, 1), size(vid, 2), length(mInds));
disp('Tracking pupil')
lastCentroid = [];

ims = zeros(size(vid,1), size(vid,2), length(mInds));
parfor i = mInds
    imgOrig = vid(:,:,i);
    img = imgOrig.*eyeMask;
    img = img < thresh & img > 0; % Threshold out the pupil
    ims(:,:,i) = imclose(img, strel('disk', 4)); % Fill in gaps
end

for i = mInds
    statusBar(i);
    imgOrig = vid(:,:,i);
    img = ims(:,:,i);
%     img = imgOrig.*eyeMask;
%     img = img < thresh & img > 0; % Threshold out the pupil
%     img = imclose(img, strel('disk', 4)); % Fill in gaps

    if usePupilArea
        CC = bwconncomp(img,4);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        if ~isempty(numPixels)
%             [biggest,idx] = max(numPixels);
            [s,idx] = sort(numPixels(:),'descend'); % sort to vector
            centroids = regionprops(CC, 'Centroid');
            centers = centroids(idx(1)).Centroid;
            radii = sqrt(s(1)/pi);
            normArea = numPixels(idx(1)) / maxArea; 
            overlayIm = imgOrig;
            overlayIm(CC.PixelIdxList{idx(1)}) = 0;
            
            %%% Make sure it region does not travel too far away. This
            %%% needs to be tested. 
%             for i = 1:length(idx)
%                 dist(i) = sqrt(sum((lastCentroid - centroids(i).Centroid).^2));
%                 overlayIm(CC.PixelIdxList{idx(i)}) = i;
%             end
%             if i == 715
%                 a = 1
%             end
            if i > 1
                if sqrt(sum((lastCentroid - centers(1,:)).^2)) > 20
                    if length(idx) > 1
                        if numPixels(idx(2))/maxArea > 0.02
                            centers = centroids(idx(2)).Centroid;
                            radii = sqrt(s(2)/pi);
                            normArea = numPixels(idx(2)) / maxArea; 
                            overlayIm = imgOrig;
                            overlayIm(CC.PixelIdxList{idx(2)}) = 0;
                        end
                    end
                end
            end
            
            

            overlayVid(:,:,i) = overlayIm;
            lastCentroid = centers(1,:);
        else
            radii = [];
            centers = [];
        end
    else
        [centers, radii, metric] = imfindcircles(img,[5 40]); 
        normArea = pi*radii(1)^2 / maxArea;
    end



    if ~isempty(radii)
        if showPupilVid
            imagesc(img); colormap gray;
            viscircles(centers(1,:), radii(1,:),'EdgeColor','r','LineWidth',0.1);
            pause(.005);
        end
        r(i) = radii(1);
        c(i,:) = centers(1,:);
        areas(i,:) = normArea(1); 
    else
        disp('cant find circle')
        if i == 1
            r = [0];
            c = [0];
        else
            r(i) = r(i-1);
            c(i,:) = c(i-1,:);
        end
    end
end


pupilA = areas;
pupilR = r;
pupilC = c;
pupilTimes = vidTimes(mInds);

%%% Generate overlay vid
if ~usePupilArea
    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',255);
    for iter = 1:length(mInds)
        i = mInds(iter);
    %     overlayVid(:,:,i) = vid(:,:,i);
        overlayVid(:,:,iter) = step(shapeInserter, vid(:,:,i), [c(iter, 1), c(iter, 2), r(iter)]);
    end
end