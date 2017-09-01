%%%% Two photon analysis script
function ImportMultiTrialDatasets2p(data)

errs = [];
for i=1:numel(data)
    try
        currData = data{i};
        mouseName = currData{1}; protocolName = currData{2};  date = currData{3}; session = currData{4};
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)];
        load(MakeBpodBehaviorPath(bpodDataPath, mouseName,  protocolName, date, session));
        
        ImportMultiTrialData2p(dataName, ...
                               MakeBpodImagePath2p(bpodImagePath, mouseName, protocolName, date, session), ...
                               SessionData)
    catch e          
      disp(e);
      errs(end+1) = i;
    end
end    


end
