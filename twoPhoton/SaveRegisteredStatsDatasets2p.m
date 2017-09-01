function SaveRegisteredStatsDatasets2p(data, varargin)
% Wrapper function to loop through different datasets. 
% Save out the max, median, variance, and skewness images across trials
% into the processed data folder. The images are taken as i.e. the maximum
% across trials of the variance within each trial. 
%
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath

p = inputParser();
p.addParameter('numTrials', [], @isnumeric); 
p.addParameter('isNLW', true, @islogical); %%% Data captured on neurolabware scope (as opposed to bruker)

p.parse(varargin{:});
numTrials = p.Results.numTrials;
isNLW = p.Results.isNLW;

errs = [];
for i=1:numel(data)
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
        SaveRegisteredStatsData2p(dataName, ...
                                 bpodTrialImagePath, ...
                                 'numTrials', numTrials, 'isNLW', isNLW, ...
                                 'zplanes', zplanes);
                                
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    