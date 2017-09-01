function SummarizeTracesDatasets2p(data, varargin)
% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for each session containing traces for each trial from the selected cells.
% Traces mat is formatted as: [cell, timepoint, trial]
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath bpodDataPath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric); 
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
tDownsampleFactor = p.Results.tDownsampleFactor;

errs = [];
for i=1:numel(data)
        close all
%     try
        currData = data{i};
        mouseName = currData{1}; protocolName = currData{2};  date = currData{3}; session = currData{4};
        if numel(currData) > 4
            zplanes = currData{5}
        else
            zplanes = [];
        end
        dataName = [mouseName '_' protocolName '_' date '_Session' num2str(session)]
        
        bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                     protocolName, date, session);
                                                 
        bpodTrialDataPath = MakeBpodBehaviorPath( bpodDataPath, mouseName, protocolName, date, session )
        
        
        nzplanes = 1;
        if ~isempty(zplanes)
            nzplanes = numel(zplanes)
        end

        for zz = 1:nzplanes
            if ~isempty(zplanes)
                zplane = zplanes(zz);
            else
                zplane = [];
            end
            
            SummarizeTracesData2p(dataName, ...
                                    bpodTrialImagePath, ...
                                    bpodTrialDataPath, ...
                                    'numTrials', numTrials, ...
                                    'isNLW', isNLW,...
                                    'zplane', zplane);
        end
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    
