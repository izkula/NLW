function PlotTrialTypes( trialTypes, success )
%PLOTTRIALTYPES(trialTypes, success)
%   vertically plot trial types
%   optional second arg color codes by whether was successful, error, or
%   3rd option

if nargin < 2
    tt = zeros(2, numel(trialTypes));
    tt(1,:) = trialTypes == 1;
    tt(2,:) = trialTypes == 2;
    imagesc(tt'); axis off; colormap hot;
else
    tt = zeros(2, numel(trialTypes));
    tt(1,trialTypes==1) = success(trialTypes == 1);
    
    tt(2,trialTypes==2) = success(trialTypes == 2);
    tt(tt==0) = -1;
    tt(1,trialTypes~=1) = 0;
    tt(2,trialTypes~=2) = 0;    
    
    imagesc(tt'); axis off; colormap(bluewhitered(256));
end

end

