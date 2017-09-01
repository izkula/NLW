function ICAsegmentationMultiTrialDatasets2p(data, varargin)
% Computes ICA components used for later cell segmentation.
%
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath

p = inputParser();
% p.addParameter('dsamp', [], @isnumeric); 

p.parse(varargin{:});
% dsamp = p.Results.dsamp;

errs = [];
for i=1:numel(data)
        close all
%     try
        currData = data{i}
        mouseName = currData{1}; protocolName = currData{2};  date = currData{3}; session = currData{4};
        if numel(currData) > 4
            zplanes = currData{5}
        else
            zplanes = [];
        end
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)]
        
        bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                     protocolName, date, session);
                                                 
        flims = ''
        nPCs = 200
        dsamp = 2
        badframes = []
        mu = 0.1
            ICAselectionMultiTrialData2p(dataName, ...
                                    bpodTrialImagePath, ...
                                    'zplanes', zplanes, ...
                                    'flims', flims, ...
                                    'nPCs', nPCs, ...
                                    'dsamp', dsamp, ...
                                    'badframes', badframes, ...
                                    'mu', mu);                                
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    

