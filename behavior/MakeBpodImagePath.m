function  output = MakeBpodImagePath( basePath, mouseName, protocolName, date, session)

if nargin < 3
    currData = mouseName;
    mouseName = currData{1}; protocolName = currData{2}; date = currData{3}; session = currData{4};
end

% output = fullfile(basePath,  'BpodImageData', mouseName, protocolName, date, ['Session' num2str(session)]);
output = fullfile(basePath, mouseName, protocolName, date, ['Session' num2str(session)]);

end

