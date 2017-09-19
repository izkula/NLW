function trials = ReshapeToTrials(alltrials, startFrames, endFrames)
%%% Reshapes from concatenated trials to a 3d matrix (cell x time x trial).

if ~exist('endFrames', 'var') || isempty(endFrames)
    endFrames = startFrames(2:end);
    startFrames = startFrames(1:end-1);
end


if size(alltrials, 2) > 1
    trials = zeros(size(alltrials, 1), ceil(max(diff(endFrames - startFrames))), numel(startFrames));
    for i = 1:numel(startFrames)
        trials(:,1:endFrames(i)-startFrames(i), i) = alltrials(:, startFrames(i):endFrames(i)-1);
    end
else
    trials = zeros(ceil(max(diff(startFrames))), numel(startFrames));
    for i = 1:numel(startFrames)
        trials(1:endFrames(i)-startFrames(i), i) = alltrials(startFrames(i):endFrames(i)-1);
    end
end


% if size(alltrials, 2) > 1
%     trials = zeros(size(alltrials, 1), ceil(max(diff(startFrames))), numel(startFrames));
%     for i = 1:numel(startFrames)-1
%         trials(:,1:startFrames(i+1)-startFrames(i), i) = alltrials(:, startFrames(i):startFrames(i+1)-1);
%     end
% else
%     trials = zeros(ceil(max(diff(startFrames))), numel(startFrames));
%     for i = 1:numel(startFrames)-1
%         trials(1:startFrames(i+1)-startFrames(i), i) = alltrials(startFrames(i):startFrames(i+1)-1);
%     end
% end



%%%%%%%ADD IN END_FRAMES?????