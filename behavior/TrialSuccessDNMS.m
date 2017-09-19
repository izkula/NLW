function successType = TrialSuccessDNMS(trialType, states)
% old function for Bpod raw data
% trialTypes tells trial types, behavior is SessionData.Trial
if isnan(states.Reward) & isnan(states.Punish)
    if trialType == 2 % should lick
        % corrrect reject
        successType = 0;
    else % shouldn't lick
        % false reject
        successType = 1;
    end
end

if ~isnan(states.Reward)
    % correct accept
    successType = 1;
end
if ~isnan(states.Punish)
    % false alarm
    successType = 0;
    trialType
end


end