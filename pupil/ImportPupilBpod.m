function  ImportPupilBpod( dataDir, varargin )
%%%% Import and process nidaq output, as well as track pupil.
%%%% dataName - i.e. 20150101/vglut1_1

%     global basePath
    
    p = inputParser();
    p.addParameter('templateDataName', '', @ischar); % If using previously computed configuration files (i.e. pupil crop coordinates, atlas, etc.)
                                                     %%% Indicate the dataName here (i.e. 20150101/vglut1_1) 
    p.addParameter('doTrackPupil', true, @islogical); % Load manta frames and track pupil
    p.addParameter('loadNidaq', true, @islogical); % Load nidaq which has time stamps of orca, manta, ball motion
    p.addParameter('doPlot', false, @islogical);
    p.addParameter('loadStimChan', false, @islogical);
    p.addParameter('loadMovement', false, @islogical);
    p.parse(varargin{:});
   
    templateDataName = p.Results.templateDataName; 
    
    doTrackPupil = p.Results.doTrackPupil;
    loadNidaq = p.Results.loadNidaq;
    doPlot = p.Results.doPlot;
    loadStimChan = p.Results.loadStimChan;
    loadMovement = p.Results.loadMovement;
    

    processedDataDir = strrep(dataDir, 'data', 'processedData');
    
    vidlist = dir(fullfile(dataDir, '*.mp4'));
    ntrials = numel(vidlist);
    
    
    allOrcaTimes = {}
    allMovementRate = {}
    allMovementTimes = {}


    allMantaTimes = {}
    allMantaFrames = {}
    allPupilA = {}
    allPupilR = {}
    allPupilC = {}
    allPupilOverlayVid = {}
    allPupilTimes = {}
    allStimOnTimes = {}
    allStimOffTimes = {}
    
    
    for i = 1:ntrials
        vidName = vidlist(i).name
        trialName = vidName(1:end-4);
        nidaqName = fullfile(dataDir, [trialName,'.bin']);
        
        
        %%%%%% Load pupil data %%%%%%%
        if doTrackPupil
            [croppedPupilFrames, thresh, eyeMask] = GetPupilFrames(dataDir, vidName, ...
                                                                'templateDataName', templateDataName, ...
                                                                'fromBpod', true);
        end

        %%%%%% Load nidaq data %%%%%%%
        if loadNidaq
            nChannels = 5; %6
            [t, ch] = LoadNidaqOutput(nidaqName, nChannels);
%             if doPlot
%                 figure, plot(t, ch); legend('1', '2', '3', '4', '5', '6');
%             end
            [mantaTimes, orcaTimes, movementRate, tMovementRate, stimTimes] = SortMantaOrcaNidaq(t, ch, doPlot, ...
                                                                                                 'doOrca', false, ...
                                                                                                 'doMovement', false, ...
                                                                                                 'doStim', true, ...
                                                                                                 'hardcodeMantaStart', false, ...
                                                                                                 'getStimOffTimes', true);
            if doTrackPupil
                mantaTimes = mantaTimes(1:min(numel(mantaTimes), size(croppedPupilFrames, 3)));
            end
            if loadStimChan
                if ~isempty(stimTimes)
                    stimOnTimes = stimTimes(1,:);
                    stimOffTimes = stimTimes(2,:);
                    fname = fullfile(processedDataDir, 'stimTimes.mat');
                    save(fname, 'stimOnTimes', 'stimOffTimes', 'orcaTimes');
                else
                    stimOnTimes = [];
                    stimOffTimes = [];
                end
            end
            disp('Done loading nidaq')
        end
        
        
        %%%%%% Do pupil tracking %%%%%%

        if doTrackPupil
            pupilTimes = mantaTimes;

            [pupilA, pupilR, pupilC, pupilTimes, pupilOverlayVid] = TrackPupil(croppedPupilFrames, ...
                                                    pupilTimes, thresh, eyeMask, 'showPupilVid', false);

            smoothFactor = 10;
            pupilA = smooth(pupilA, smoothFactor);
            pupilR = smooth(pupilR, smoothFactor); 
            if doPlot
                figure, subplot(1,2,1); plot(pupilTimes, pupilR), title('radius')
                subplot(1,2,2); plot(pupilTimes, sqrt(sum(pupilC.^2, 2))), title(['pupil c: ', dataName])
                drawnow
                pause(0.1);
            end
        %     figure, VideoSlider(pupilOverlayVid)
        else
            croppedPupilFrames = [];
            mantaTimes = [];
        end

        
        allOrcaTimes{i} = orcaTimes
        allMovementRate{i} = movementRate;
        allMovementTimes{i} = tMovementRate;

        allMantaTimes{i} = mantaTimes';
        allPupilA{i} = pupilA;
        allPupilR{i} = pupilR;
        allPupilC{i} = pupilC;
        allPupilTimes{i} = pupilTimes';
        allStimOnTimes{i} = stimOnTimes;
        allStimOffTimes{i} = stimOffTimes;
        
        allPupilOverlayVid{i} = pupilOverlayVid;
        allMantaFrames{i} = croppedPupilFrames;
    end


%%%%% Save out covariates %%%%%%
%     fname = fullfile(processedDataDir, 'covariates.h5');
%     if exist(fname, 'file')
%         disp(['deleting old ' fname]); delete(fname);
%     end
%     
%     if ~isempty(orcaTimes)
%         H5MakeNewDataset(fname, '/orcaTimes', cell2mat(allOrcaTimes));
%     end
%     
%     if ~isempty(movementRate)
%         H5MakeNewDataset(fname, '/movementRate',
%         cell2mat(allMovementRate));
%         H5MakeNewDataset(fname, '/movementTimes', cell2mat(allMovementTimes));
%     end
% 
%     if doTrackPupil        
%         H5MakeNewDataset(fname, '/mantaTimes', cell2mat(allMantaTimes));        
%         H5MakeNewDataset(fname, '/pupilA', cell2mat(allPupilA));
%         H5MakeNewDataset(fname, '/pupilR', cell2mat(allPupilR));
%         H5MakeNewDataset(fname, '/pupilC', cell2mat(allPupilC));
%         H5MakeNewDataset(fname, '/pupilTimes', cell2mat(allPupilTimes));   %%% Note: times should be relative to beginning of trial when Nidaq starts recording
%                                                                  %%%   (i.e. mantaTimes do not start at 0)
% %         H5MakeNewDataset(fname, '/mantaFrames', allMantaFrames );
% %         H5MakeNewDataset(fname, '/pupilOverlayVid', allPupilOverlayVid);
%     end
% 
%     if loadStimChan
%         H5MakeNewDataset(fname, '/stimOnTimes', cell2mat(allStimOnTimes));        
%         H5MakeNewDataset(fname, '/stimOffTimes', cell2mat(allStimOffTimes));        
%     end

    save(fullfile(processedDataDir, 'pupilTraces.mat'), 'allMantaTimes', 'allPupilA', 'allPupilR', 'allPupilC', 'allPupilTimes', 'allStimOnTimes', 'allStimOffTimes');
    
    save(fullfile(processedDataDir, 'pupilVids.mat'), 'allPupilOverlayVid', '-v7.3');
    end




