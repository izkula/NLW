function  output = MakeBpodImagePath2p( basePath, mouseName, protocolName, date, session)

output = fullfile(basePath,  mouseName, protocolName, date, ['Session' num2str(session)]);

end