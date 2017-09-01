function SaveTracesMultiTrialDatasets2p(data, varargin)
% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for each session containing traces for each trial from the selected cells.
% Traces mat is formatted as: [cell, timepoint, trial]
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric); 
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('isContinuousTrial',false,@islogical);
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
isContinuousTrial = p.Results.isContinuousTrial;
tDownsampleFactor = p.Results.tDownsampleFactor;

errs = [];
for i=1:numel(data)
        close all
%     try
        currData = data{i};
        mouseName = currData{1}; protocolName = currData{2};  date = currData{3}; session = currData{4};
        zplanes = currData{5};
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)]
        
        bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                     protocolName, date, session);
                                                 
        if isContinuousTrial
            SaveTracesContMultiTrialData2p(dataName, ...
                                        bpodTrialImagePath, ...
                                        'numTrials', numTrials, ...
                                        'isNLW', isNLW, ...
                                        'zplanes', zplanes, ...
                                        'tDownsampleFactor', tDownsampleFactor);
        else
            SaveTracesMultiTrialData2p(dataName, ...
                                        bpodTrialImagePath, ...
                                        'numTrials', numTrials, ...
                                        'isNLW', isNLW, ...
                                        'zplanes', zplanes);
        end
                                
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    
