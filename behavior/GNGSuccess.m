function successType = GNGSuccess(trialType, behavior)
% old function for Bpod raw data
% trialTypes tells trial types, behavior is SessionData.Trial
if isnan(behavior.States.Reward) & isnan(behavior.States.Punish)
    if trialType == 1 % should lick
        % corrrect reject
        successType = 0;
    else % shouldn't lick
        % false reject
        successType = 1;
    end
end

if ~isnan(behavior.States.Reward)
    % correct accept
    successType = 1;
end
if ~isnan(behavior.States.Punish)
    % false alarm
    successType = 0;
end


end