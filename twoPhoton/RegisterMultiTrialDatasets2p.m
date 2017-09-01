function RegisterMultiTrialDatasets2p(data, varargin)
% Motion correct 2p datasets saved out using bpod
%
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2, []}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1, []}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum, [zplanes]}}
%

global bpodImagePath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric);
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)
p.addParameter('useTemplate', true, @islogical); %%% Register all trials to the first trial
p.addParameter('doDebug', false, @islogical); %%% Register all trials to the first trial
p.addParameter('isContinuousTrial',false,@islogical);
p.addParameter('tDownsampleFactor', 3, @isnumeric);

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;
useTemplate = p.Results.useTemplate;
doDebug = p.Results.doDebug;
isContinuousTrial = p.Results.isContinuousTrial;
tDownsampleFactor = p.Results.tDownsampleFactor;

errs = [];
allStartTime = tic;
for i=1:numel(data)
     try
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
        if ~isNLW                                     
            FixTrialNames(bpodTrialImagePath);
        end
        

        RegisterMultiTrialData2p(dataName, ...
                                 bpodTrialImagePath, ...
                                 'numTrials', numTrials, 'isNLW', isNLW, ...
                                 'useTemplate', useTemplate, 'zplanes', zplanes, ...
                                 'doDebug', doDebug, ...
                                 'isContinuousTrial', isContinuousTrial, ...
                                  'tDownsampleFactor', tDownsampleFactor);
        
        disp('*******************************************************')
        disp('*******************************************************')
        disp(['Time of finishing ', dataName])    
        disp(toc(allStartTime))
        disp('*******************************************************')
        disp('*******************************************************')
     catch e          
       disp(e);
       errs(end+1) = i;
     end
     disp(errs);
end    




