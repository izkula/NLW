function ICAsegmentationMultiTrialDatasets2p(data, varargin)
% Cell segmentation based on ICA components. 
%
% Parameters
% data: a cell array, i.e. {{'vglut1m15', 'Image2PShapeOlfGNG', 'Dec10_2015', 2}, ...
%                           {'vglut1m16', 'Image2PShapeOlfGNG', 'Dec10_2015', 1}}, ...
%                           {'mouseName', 'ProtocolName', 'Mon00_year',  sessNum}}
%

global bpodImagePath

p = inputParser();
p.addParameter('doPlotting', 0, @isnumeric); 

p.parse(varargin{:});
plotting = p.Results.doPlotting;

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
                                                 
        smwidth = 0.2
        thresh = 4
        arealims = 20
    
        ICAsegmentationMultiTrialData2p(dataName, ...
                                    bpodTrialImagePath, ...
                                    'smwidth', smwidth, ...
                                    'thresh', thresh, ...
                                    'arealims', arealims, ...
                                    'plotting', plotting, ...
                                    'zplanes', zplanes);                                
%     catch e          
%       disp(e);
%       errs(end+1) = i;
%     end
end    

