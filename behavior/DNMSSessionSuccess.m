function success = DNMSSessionSuccess( SessionData, isStruct)
%GNGSESSIONSUCCESS compute successes and failures on each trial from
%SessionData

if nargin < 2
    isStruct = true;
end
success = [];
for i=1:SessionData.nTrials
    if isStruct
        success = [success TrialSuccessDNMS(SessionData.TrialTypes(i), SessionData.states{i})];
    else
%         success = [success GNGSuccess(SessionData.TrialTypes(i), SessionData.RawEvents.Trial{i})];
    end
end


end

