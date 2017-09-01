function SelectCellsICAMultiTrialDatasets2p(data, varargin)
% Provides a gui interface for selecting cells for each of the specified
% datasets (each of which consists of multiple trials).
%
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric); %%% Irrelevant for manual selection, since we currently only use the first trial.
p.addParameter('automateSelection', true, @islogical); %%% Use automated cell selection followed by manual adjustment
p.addParameter('cellRadius', [], @isnumeric);
p.addParameter('doAppendTiffs', false, @islogical);
p.addParameter('isNLW', true, @islogical);


p.parse(varargin{:});
numTrials = p.Results.numTrials;
automateSelection = p.Results.automateSelection;
cellRadius = p.Results.cellRadius;
doAppendTiffs = p.Results.doAppendTiffs;
isNLW = p.Results.isNLW;

errs = [];
for i=1:numel(data)
        close all
%     try
        currData = data{i}
        mouseName = currData{1}; protocolName = currData{2};  date = currData{3}; session = currData{4};
        zplanes = currData{5};
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)]
        
        bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                     protocolName, date, session);
        SelectCellsICAMultiTrialData2p(dataName, ...
                                    bpodTrialImagePath, ...
                                    'numTrials', numTrials, ...
                                    'automateSelection', automateSelection, ...
                                    'cellRadius', cellRadius, ...
                                    'zplanes', zplanes, ...
                                    'doAppendTiffs', doAppendTiffs, ...
                                    'isNLW', isNLW);
                                
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    

