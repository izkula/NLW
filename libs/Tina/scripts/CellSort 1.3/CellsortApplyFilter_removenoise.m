function cell_sig = CellsortApplyFilter_removenoise(fn, ica_segments, flims, subtractmean)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters
%
% Inputs:
%     fn - file name of TIFF movie file
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractmean - boolean specifying whether or not to subtract the mean
%     fluorescence of each time frame
%
% Outputs:
%     cell_sig - cellular signals
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic

if (nargin<2)||isempty(flims)
    nt = tiff_frames(fn);
    flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<4
    subtractmean = 1;
end

[pixw,pixh] = size(imread(fn,1));
k=0;

cell_sig = zeros(size(ica_segments,1), nt);
ica_segments = reshape(ica_segments, [], pixw*pixh);


fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
while k<nt
    ntcurr = min(500, nt-k);
    mov = zeros(pixw, pixh, ntcurr);
    for j=1:ntcurr
        movcurr = imread(fn, j+k+flims(1)-1);
        mov(:,:,j) = movcurr;
    end
    %mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean

    if subtractmean
        % Subtract the mean of each frame
         movtm = mean(mean(mov,1),2);
         mov = mov - repmat(movtm,[pixw,pixh,1]);
    end

    mov = reshape(mov, pixw*pixh, ntcurr);
    cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;

    k=k+ntcurr;
    fprintf('Loaded %3.0f frames; ', k)
    toc
end

% fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
% while k<nt
%     ntcurr = min(500, nt-k);
%     mov = zeros(pixw, pixh, ntcurr);
%     for j=1:ntcurr
%         movcurr = imread(fn, j+k+flims(1)-1);
%         mov(:,:,j) = movcurr;
%     end
%     %mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean
%     
%     if subtractmean
%         % Subtract the mean of each frame
%         %         figure(100);
%         %         plot(movtm);
%         %         thresh=input('enter threshold for flickering');
%         %         hold on
%         %         plot([0 nt],[thresh thresh],'-r');
%         %         mov (:,:,movtm>thresh) = bsxfun(@minus,mov(:,:,movtm>thresh),mean(mov(:,:,movtm>thresh),3));
%         %         mov (:,:,movtm<thresh) = bsxfun(@minus,mov(:,:,movtm<thresh),mean(mov(:,:,movtm<thresh),3));
%         movtm = mean(mean(mov,1),2);
%         mov = mov - repmat(movtm,[pixw,pixh,1]);
%     end
% 
%     mov = reshape(mov, pixw*pixh, ntcurr);
%     movm = mean(mov,2); % Average over time
%     movm = reshape(movm, pixw, pixh);
%     cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;
%     %movtm(:,k+[1:ntcurr])=mean(mov,1);
% 
%     k=k+ntcurr;
%     fprintf('Loaded %3.0f frames; ', k)
%     toc
%     
% end

    function j = tiff_frames(fn)
        %
        % n = tiff_frames(filename)
        %
        % Returns the number of slices in a TIFF stack.
        %
        % Modified April 9, 2013 for compatibility with MATLAB 2012b

        j = length(imfinfo(fn));
