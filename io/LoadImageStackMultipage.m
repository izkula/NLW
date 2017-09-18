function currImages = LoadImageStackMultipage( fname, subSeq)
% Load video stack from a multipage tiff.


info = imfinfo(fname);
numImages = numel(info);
im =  imread(fname, 1, 'Info', info);
imSize = size(im);
if ~exist('subSeq', 'var') || isempty(subSeq)
    subSeq = 1:numImages;
end

currImages = zeros(imSize(1), imSize(2), length(subSeq));
iter = 1;
for k = subSeq
    try
        currImages(:,:,iter) = imread(fname, k, 'Info', info);
    catch e
    end
    iter = iter + 1;
end
