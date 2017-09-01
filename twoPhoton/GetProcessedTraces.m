function C = GetProcessedTraces(traces, varargin)
%%% For dex 2p data

p = inputParser();
p.addParameter('doPlot',false, @islogical);
p.addParameter('doDeconvolve', false, @islogical)

p.parse(varargin{:});
doPlot = p.Results.doPlot;
doDeconvolve = p.Results.doDeconvolve;

%% Remove NaN traces and concatenate trials
B = reshape(traces, [size(traces,1), size(traces,2)*size(traces,3)]); %%% Append all minute vids
labels = 1:size(B, 1);
labels(~any(~isnan(B), 2)) = [];
B(~any(~isnan(B), 2),:)=[];
B(var(B, 0, 2)>.1,:)=[]; %%% Only use traces with a minimum variance (i.e. amount of activity)
B(min(B, [], 2)<-.6,:)=[];

if doPlot
    figure, imagesc(B)
end
% figure, plot(var(B, 0,2))



%% Deconvolve traces
% Says will: not sure what most of these really mean...

if doDeconvolve
    C = zeros(size(B));
    for kk = 1:size(B, 1)
        kk
    options = struct();
    options.dt = 1/30.;
    options.MaxIter = 10.;
    options.MaxInerIter = 50;
    options.TauStd = [0.2, 2];
    options.default_g = [0.6,0.9];
    options.method = 'cvx';
    b = [];
    c1 = [];
    g = [];
    sn = [];

    trace = B(kk, :);
    [c,~,~,~,~,sp] = MCEM_foopsi(trace,b,c1,g,sn,options);

    B(kk,:) = B;
    end
end
%% Find spike times
C = B; %%% C will contain processed spikes

threshMult = 3;
for i = 1:size(B,1)
    cellTrace = C(i,:);
    sCellTrace = smooth(cellTrace, 10, 'sgolay');
%     figure, plot(cellTrace, 'r'); hold on; plot(sCellTrace, 'b')
%     scaleFactor = robustfit(FF, smooth(cellTrace', 5));
% %      plot(smooth(cellTrace, 5), 'r'); hold on; plot(scaleFactor(1) + FF*scaleFactor(2), 'b')
%     cellTraceNew = smooth(cellTrace, 5) - (scaleFactor(1) + FF*scaleFactor(2));
%     plot(smooth(cellTrace,5), 'r'); hold on; plot(cellTraceNew, 'b'); hold off;
%      pause(1)
    sCellTrace(sCellTrace < threshMult*std(sCellTrace)) = 0;
    C(i,:) = sCellTrace;
end

if doPlot
    figure, imagesc(B(:, 1:300))
    figure, imagesc(C(:, 1:300))
end

end