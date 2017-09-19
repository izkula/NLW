function s = Bpod2Struct( SessionData )
% convert Bpod data to struct with events and times
s = struct();
s.TrialStartTimeStamp = SessionData.TrialStartTimestamp;
s.nTrials = SessionData.nTrials;
if isfield(SessionData, 'TrialTypes')
    s.TrialTypes = SessionData.TrialTypes;
end
if isfield(SessionData, 'MatchTypes')
    s.MatchTypes = SessionData.MatchTypes;
end
s.events = cell(s.nTrials,1);
s.states = cell(s.nTrials,1);
for i=1:SessionData.nTrials
    currTrial = SessionData.RawEvents.Trial{i};
    s.events{i} = currTrial.Events;
    s.states{i} = currTrial.States;    
end
end

