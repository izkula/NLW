function [ rxnTimes ] = ReactionTime( eventTimes, timeInterval )
%REACTIONTIME 
% Input: eventTimes (cell array of licks), timeInterval ([t1,t2) in seconds)
% Get the time of the first event within timeInterval 
rxnTimes = zeros(numel(eventTimes),1);
for i=1:numel(eventTimes)
    currEvents = eventTimes{i};
    if ~isempty(currEvents)
        currSubset = currEvents(currEvents >= timeInterval(1) & currEvents < timeInterval(2));
        if ~isempty(currSubset) 
            rxnTimes(i) = min(currSubset);
        end
    end
end

end

