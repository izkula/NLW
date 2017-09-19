function eventTimes = GetEventTimes( SessionData, eventName, isStruct )
%GETEVENTTIMES Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    isStruct = true;
end

eventTimes = {};
for i=1:SessionData.nTrials
    if isStruct
        currTrial = SessionData.events{i};
    else
        currTrial = SessionData.RawEvents.Trial{i}.Events;
    end
    if isfield(currTrial, eventName)
        eventTimes{i} = getfield(currTrial, eventName);
    else
        eventTimes{i} = [];
    end
end

end

