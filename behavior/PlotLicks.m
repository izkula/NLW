function PlotLicks( SessionData, showSuccess, isStruct)
if nargin < 3
    isStruct = true;
end
if showSuccess
    success = GNGSessionSuccess(SessionData);
end
for i=1:SessionData.nTrials
    if isStruct
        currEvents = SessionData.events{i};
    else
        currEvents = SessionData.RawEvents.Trial{i}.Events;
    end
    
    if isfield(currEvents, 'Port1In')                    
        if showSuccess
            if success(i) == 1
                plot(currEvents.Port1In, i, '.k'); hold on;
            else
                plot(currEvents.Port1In, i, '.r'); hold on;
            end
        else
            plot(currEvents.Port1In, i, '.k'); hold on;
        end
        if isStruct
            if ~any(isnan(SessionData.states{i}.Reward))
                plot(SessionData.states{i}.Reward(1), i, 'go'); hold on;
            end
        else
            if ~any(isnan(SessionData.RawEvents.Trial{i}.States.Reward))
                plot(SessionData.RawEvents.Trial{i}.States.Reward(1), i, 'go'); hold on;
            end
        end
            
    end
end
ylim([0 SessionData.nTrials]);
end

