function [mantaTimes, orcaTimes, movementRate, tMovementRate, stimTimes] = SortMantaOrcaNidaq(t, ch, doPlot, varargin)
%%% Recovers the manta and orca frame times based on the Nidaq timeseries
%%% (as loaded using LoadNidaqOutput). Also optionally outputs a time
%%% series of the mouse's movement, based on the optical encoder tracking
%%% the running ball. 
%%% Make sure that the channel settings are correct at the beginning of
%%% this file. 

if nargin < 3
    doPlot = false;
end

p = inputParser();
p.addParameter('doMovement', true, @islogical); % If using previously computed configuration files (i.e. pupil crop coordinates, atlas, etc.)
                                                 %%% Indicate the dataName here (i.e. 20150101/vglut1_1) 
p.addParameter('doOrca', true, @islogical); % Load manta frames and track pupil
p.addParameter('doStim', false, @islogical);
p.addParameter('hardcodeMantaStart', true, @islogical); %%% To deal with a bug
p.addParameter('getStimOffTimes', false, @islogical); %%% To deal with a bug
p.parse(varargin{:});

doMovement = p.Results.doMovement;
doOrca = p.Results.doOrca;
doStim = p.Results.doStim;
hardcodeMantaStart = p.Results.hardcodeMantaStart;
getStimOffTimes = p.Results.getStimOffTimes;

movementChan = 4;
orcaChan = 6;
mantaChan = 1;
stimChan = 2;

% figure, plot(t, ch(movementChan, :)), title('movement');
% figure, plot(t, ch(orcaChan, :)), title('orca');
% figure, plot(t, smooth(ch(mantaChan, :), 3)), title('manta');

doUseSmooth = false;

if doUseSmooth
    aa = ch(mantaChan, :);
    aa = smooth(aa, 100);
    aa = double(aa>mean(aa));
    bb = find(aa(100:end)>0);
    ind = bb(1);
else
%%% Find first manta onset
    [~, ind] = max(diff(ch(mantaChan, :)));
end
    
if hardcodeMantaStart
    ind = 38400; %%%% This line is a temporary fix when manta frames weren't recorded. Not sure what it impacts in previous analysis. 
end

%%% Get manta frame times
% figure()
if doUseSmooth
    [mantaPks, mantaTimes] = findpeaks(aa(ind:end),t(ind:end),'MinPeakDistance',.01);
else
    [mantaPks, mantaTimes] = findpeaks(ch(mantaChan, ind:end),t(ind:end),'MinPeakDistance',.01);
end

%%% Get orca frame times
if doOrca
    [orcaPks, orcaTimes] = findpeaks(diff(ch(orcaChan, ind:end)),t(ind:end-1),'MinPeakDistance',.01, 'MinPeakHeight', 1);
    if orcaTimes(2) - orcaTimes(1) > 0.5
        orcaPks = orcaPks(2:end-1);
        orcaTimes = orcaTimes(2:end-1);
    end
else
    orcaPks = [];
    orcaTimes = [];
end

%%% Get stim times
if doStim
    stimTrace = ch(stimChan, ind:end);
    stimTrace = SharpSmooth(abs(stimTrace), .5*max(stimTrace), 10);
    [stimOnPks, stimOnTimes] = findpeaks(diff(stimTrace),t(ind:end-1),'MinPeakDistance',.01, 'MinPeakHeight', .5*max(stimTrace(:)));
    [stimOffPks, stimOffTimes] = findpeaks(-diff(stimTrace),t(ind:end-1),'MinPeakDistance',.01, 'MinPeakHeight', .5*max(stimTrace(:)));
    if numel(stimOnTimes) > numel(stimOffTimes)
        stimOnTimes(end) = [];
    end
    if numel(stimOffTimes) > numel(stimOnTimes)
        stimOffTimes(1) = [];
    end
%     stimTimes = [ stimOnTimes; stimOffTimes];
    if getStimOffTimes
        stimTimes = [stimOnTimes; stimOffTimes]; %%%% THIS IS A HACK TO FIX A BUG!!!
    else
        stimTimes = [stimOnTimes; stimOnTimes]; %%%% THIS IS A HACK TO FIX A BUG!!!
    end
else
    stimTimes = [];
end

if doMovement
    %%% Get movement rate
    %%% Rotary encoder makes 1024 pulses per revolution (using the E6C2-CWZ3E). 
    % The ball is 6 inches in diameter. Thus, each pulse corresponds to the
    % mouse walking 0.467 mm.

    % velocity = VelocityFromRotaryEncoder(ch(movementChan, ind:end), mantaTimes(2) - mantaTimes(1), 6, 1024);
    velocity = VelocityFromRotaryEncoder(ch(movementChan, ind:end), t, mantaTimes(2) - mantaTimes(1), 6, 1024);
    warning('Need to test VelocityFromRotaryEncoder in SortMantaOrcaNidaq')

    %%% Original. Refactor has not been tested. 
    % distPerPulse = 6*25.4*pi/(1024*2); % 6 inch diameter ball. 1024 pulses per revolution, and an up and down voltage for each step. 
    % pulses = abs(diff(ch(movementChan, ind:end)));
    % pulses(pulses < 0.1) = 0; 
    % pulses(pulses > 0) = 1;
    % position = cumsum(pulses)*distPerPulse;
    % dtManta = mantaTimes(2) - mantaTimes(1); 
    % dtEncoder = t(2) - t(1); 
    % indPerManta = round(dtManta/dtEncoder); % Should be an exact multiple of one another
    % velocity = zeros(size(mantaTimes));
    % p = position(1:indPerManta:end);
    % velocity = (p(2:end) - p(1:end-1))/dtManta;

    if numel(velocity) > numel(mantaTimes)
        velocity = velocity(1:numel(mantaTimes));
    elseif numel(velocity) < numel(mantaTimes)
        velocity = [zeros(1, numel(mantaTimes) - numel(velocity)), velocity];
    end

    movementRate = velocity; 
    tMovementRate = mantaTimes; 
    % smoothWin = 501;
    % movementRate = smooth(abs(diff(ch(movementChan, ind:end), smoothWin)));
    % movementRate = 10*movementRate/max(movementRate(:));
    % tMovementRate = t(ind+(smoothWin-1)/2:(smoothWin-1)/2 + ind+length(movementRate)-1);
    % 
    % 
    % decimationFactor = round(length(movementRate)/length(mantaTimes)/2)
    % movementRate = decimate(movementRate, decimationFactor);
    % tMovementRate = decimate(tMovementRate, decimationFactor);
    % movementRate(movementRate < 0.1) = 0;
else
    movementRate = [];
    tMovementRate = [];
end

%%% Plot all waveforms
if doPlot
    figure(), plot(mantaTimes, mantaPks, 'go');
    hold on; plot(t(ind:end), ch(mantaChan, ind:end))
    hold on, plot(t(ind), 5, 'ro')
    hold on; plot(t(ind:end), ch(orcaChan, ind:end), 'k-')
    if doOrca
        hold on; plot(orcaTimes, orcaPks, 'mo')
    end
    if doMovement
        hold on; plot(t(ind:end-1), abs(diff(ch(movementChan, ind:end))))
        hold on; plot(tMovementRate, movementRate)
    end
end