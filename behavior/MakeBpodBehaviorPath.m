function output = MakeBpodBehaviorPath( basePath, mouseName, protocolName, date, session )
fname = [mouseName '_' protocolName '_' date '_' 'Session' num2str(session) '.mat'];
output = fullfile(basePath, mouseName, protocolName, 'Session Data', fname);

end

