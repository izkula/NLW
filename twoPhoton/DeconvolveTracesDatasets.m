function DeconvolveTracesDatasets( data, varargin )

% After SelectCellsMultiTrialDatasets has been called, this function
% saves out a mat for each session containing traces for each trial from the selected cells.
% Traces mat is formatted as: [cell, timepoint, trial]
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global basePath bpodImagePath resultsPath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric); 
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('tDownsampleFactor', [], @isnumeric); 

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
tDownsampleFactor = p.Results.tDownsampleFactor;

errs = [];
for i=1:numel(data)
%     try
        currData = data{i};

        mouseName = currData{1}; protocolName = currData{2};  currDate = currData{3}; session = currData{4}; 
        if numel(currData) > 4
            zplanes = currData{5}
        else
            zplanes = [];
        end
        
        nzplanes = 1;
        if ~isempty(zplanes)
            nzplanes = numel(zplanes)
        end
        
        tDownsampleFactor = tDownsampleFactor/nzplanes;
        
        dt = 1/30.;
        dt = dt*tDownsampleFactor;

        for zz = 1:nzplanes
            if ~isempty(zplanes)
                zplane = zplanes(zz);
            else
                zplane = [];
            end
            
            maxT = [];
            dataset = Dataset2P(mouseName, protocolName, currDate, session, dt, maxT, isNLW, zplane); %%% A big struct that you save out
            bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                         protocolName, currDate, session);        

            saveDir = strrep(bpodTrialImagePath, bpodImagePath, fullfile(resultsPath, '2p/'))
            mkdir(saveDir);
            if ~isempty(zplane)
                zstr = ['_z', num2str(zplane), '_'];
            else
                zstr = '';
            end
            save(fullfile(saveDir, ['deconvolved', zstr, '.mat']), 'dataset');
        end

%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    

