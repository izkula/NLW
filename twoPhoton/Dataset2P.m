classdef Dataset2P
    %DATASET2P Two photon dataset 
    %  Includes traces, image, and behavior
    
    properties
        Ncells = 0;
        Nt = 0;
        Ntrials = 0;
        img = [];
        behavior = [];
        mouseName = '';
        protocolName = '';
        saveDir = '';
        date = '';
        session = -1;
        success = [];
        dt = 1/30.;
        goodTraces = [];
        isNLW = false;
        zplane = [];
        metaDatas = [];
        useTrial = []; %%% Indicates whether bpod properly recorded start frame time
    end
    
    methods
        function obj = Dataset2P(mouseName, protocolName, date, session, dt, maxT, isNLW, zplane)
            global bpodImagePath bpodDataPath resultsPath
           
            bpodTrialImagePath =  MakeBpodImagePath2p(bpodImagePath, mouseName, ...
                                                     protocolName, date, session);                                                 
            bpodTrialDataPath = MakeBpodBehaviorPath( bpodDataPath, mouseName, protocolName, date, session );
            obj.mouseName = mouseName;
            obj.protocolName = protocolName;
            obj.date = date;
            obj.session = session;
            obj.dt = dt;
            
            %%% Load traces into matlab                               
            obj.saveDir = fullfile(resultsPath, '2p', mouseName, protocolName, date, ['Session' num2str(session)]); %% i.e. ~/Dropbox/oeg_results/2p/gad2m11/Image2POlfGoNoGo/Dec16_2015/Session1
            if exist('isNLW', 'var') || ~isempty(isNLW)
                obj.isNLW = isNLW;
            end
            if exist('zplane', 'var') || ~isempty(zplane)
                obj.zplane = zplane;
                zstr = ['_z', num2str(zplane), '_'];
            else
                obj.zplane = [];
                zstr = '';
            end
            
            %%% Load traces
            S = load(fullfile(obj.saveDir, ['traces', zstr, '.mat']));
            obj.img = S;
            
            %%% If NLW, adjust framerate if multiple planes
            if isNLW
                obj.useTrial = S.useTrial;
                obj.metaDatas = S.metaDatas;
                try
                    if S.metaDatas{1}.volscan
                        nZPlanesRecorded = numel(S.metaDatas{1}.otwave_um);
                    else
                        nZPlanesRecorded = 1;
                    end
                catch
                    nZPlanesRecorded = 1;
                end
                disp(['nZPlanesRecord: ', num2str(nZPlanesRecorded)]);
                dt = dt*nZPlanesRecorded; 
            end

            %%% Load bpod file
            copyfile(bpodTrialDataPath,fullfile(obj.saveDir, 'behavior.mat'));
            D = load(bpodTrialDataPath);
            s = Bpod2Struct(D.SessionData);          
            obj.behavior = s;
            [Ncell, Ntime, Ntrial] = size(S.traces);      
            
            %%% Determine frame times
            try
                frameTimes = s.events{1}.BNC2High; % Just use timestamps from the first trial. 
            catch
                t_per_frame = dt;
                warning(['No frame times found; using t_per_frame = ', num2str(t_per_frame)])
                frameTimes = [0:size(S.fullField, 1)]*t_per_frame; 
            end
            if exist('maxT', 'var') || ~isempty(maxT)
               maxT =  min(size(S.traces,2), numel(frameTimes)); % Number of 2p frames
            end
            frameTimes = frameTimes(1:maxT); % We double checked that the 2p outputs *extra* frame TTLs at the *end*.
            obj.img.frameTimes = frameTimes;
            

            %obj.Nt = size(S.traces, 2); % Number of 2p frames
            obj.Nt = maxT; % Number of 2p frames
            obj.Ntrials = size(S.traces, 3);

            %%% Load neuropil traces
            if iscell(S.neuropilTraces)
                neuropil = zeros(size(S.neuropilTraces{1},1),obj.Nt,numel(S.neuropilTraces));
                for i=1:numel(S.neuropilTraces)
                    NN = S.neuropilTraces{i};
                    if size(NN, 2) < obj.Nt
                        break
                    end
                    neuropil(:,:,i) = NN(:,1:obj.Nt);
                end
            else
                neuropil = S.neuropilTraces;
            end
            obj.img.neuropilTraces = neuropil;

            
            %%%% Potentially remove this legacy code
            %df = obj.ComputeDF(obj.img.neuropilTraces);
            %dff = obj.ComputeDFF(obj.img.traces-(obj.img.neuropilTraces-df), 30);
            %traces = reshape(dff, [Ncell, Ntime-30+1, Ntrial]);            
            %dfNeuropil = zeros(size(obj.img.neuropilTraces));
            %for i=1:size(dfNeuropil,3)
                %dfNeuropil(:,:,i) = obj.img.neuropilTraces(:,:,i) ./ repmat(squeeze(median(obj.img.neuropilTraces(:,:,i),2)),1,size(obj.img.neuropilTraces,2))-1.0;                
            %end
            
            %%% Remove dim cells
            goodTraces = [];
            for i=1:size(obj.img.traces,1)
                currTrace = obj.img.traces(i,:,:);
                currNeuropil = obj.img.neuropilTraces(i,:,:);
                if mean(currTrace(:)) >= 1.05*mean(currNeuropil(:))
                    goodTraces = [goodTraces i];
                end
            end
            obj.goodTraces = goodTraces;
            disp('Good traces: ')
            disp(goodTraces)
            obj.Ncells = numel(goodTraces);
           
            
            %%% Subtract neuropil from neural trace
            alpha = 0.7;
            dff = obj.ComputeRunningDFF(obj.img.traces(obj.goodTraces,1:maxT,1:obj.Ntrials));
            neuropilDff = obj.ComputeRunningDFF(obj.img.neuropilTraces(obj.goodTraces,1:maxT,1:obj.Ntrials));
            dff = dff - alpha*neuropilDff;
            [dff, deconvCa, spikes] = obj.Deconvolve(dff, dt, 1*round(1/dt)); % ignore first 1 s
            obj.img.traces = dff;
            obj.img.deconv = deconvCa;
            obj.img.spikes = spikes;
            %%%c = 2; figure, plot(deconvCa(c,:,1)'); hold on;plot(smooth(dff(c,:,1)')); hold on; plot(spikes(c, :, 1)'); %% For plotting
            for i=1:size(obj.img.fullField,2)
                temp = obj.img.fullField(:,i);
                temp(isnan(temp)) = 0;
                obj.img.fullField(:,i) = temp; %./ repmat(median(temp,1), size(temp,1),1) - 1.0;
            end
            %obj.dffTraces = obj.ComputeDFF(30); % remove first second
        end
               
        function df = ComputeDF(obj, t)
            [nCell, nTotalTime]  = size(t);

            smoothedTSeries = t;
            df = zeros(size(smoothedTSeries));
            
            for i=1:nCell
                smoothedTSeries(i,:) = smooth(smoothedTSeries(i,:), 30, 'lowess');
                df(i,:) = smooth(smoothedTSeries(i,:), 1800, 'moving');
            end
        end
        
        function dff = ComputeDFF(obj, t, startPos)
            if nargin < 2
                startPos = 1;
            end
            perCellTSeries = obj.GetConcatenatedTimeseries(t, startPos);
            [nCell, nTotalTime]  = size(perCellTSeries);
            df = obj.ComputeDF(perCellTSeries);
            dff = perCellTSeries ./ df - 1.0;
        end
        
        function output = ComputeRunningDFF(obj, traces)
            % traces in [ncells, ntime, ntrials]
            [Ncell, Ntime, Ntrial] = size(traces);
            output = zeros(Ncell, Ntime, Ntrial);
            for i=1:size(traces,1)
                for j=1:size(traces,3)
                    temp = squeeze(traces(i,:,j));
                    %temp = reshape(mean(reshape(temp,2,[]),1),[],size(temp,2))'; % downsample each trace 2x

%                     temp = temp' ./ running_percentile(temp, numel(temp), 20) - 1.0;
                    baseline = prctile(temp(temp>0), 20);
                    temp = temp' ./ baseline - 1.0;
                    

                    output(i,:,j) = temp;
                    
                    %traces(i,:,j) = temp' ./ smooth(temp, 1800, 'moving') - 1.0;
                    %traces(i,:,j) = temp' ./ median(temp) - 1.0;
                end
            end
        end
        
        function output = GetZScoreData(obj)
            output = zeros(size(obj.img.deconv));
            [Ncell, Ntime, Ntrial] = size(obj.img.deconv);
            for i=1:Ncell
                for j=1:Ntrial
                    output(i,:,j) = zscore(squeeze(obj.img.deconv(i,:,j)));
                end
            end
        end
        
        function [go, nogo] = GetTrialTypeAvg(obj,doZscore,justActive)
            
            if doZscore
                traces = obj.GetZScoreData();
            else
                traces = obj.img.deconv;
            end
            if justActive
                traces = traces(obj.img.taskLabels,:,:);
            end
            go = mean(traces(:,:,obj.behavior.TrialTypes==1 & obj.success == 1),3);
            nogo = mean(traces(:,:,obj.behavior.TrialTypes==2 & obj.success == 1),3);
            
        end
        
        function [dists, cc] = AllPairsDistsAndCC(obj)
            dists = pdist(obj.img.cellCentroids(obj.goodTraces,:), 'euclidean')';
            flatData = reshape(obj.img.traces, [obj.Ncells, obj.Nt*obj.Ntrials]);
            cc = corr(flatData');
            cc = GetLowerTriVals(cc);
        end
        
        function dff = GetDffFullFieldTraces(obj)
            % get zscored, dff full field traces downsampled 30 hz -> 15 hz and DFF
            [T N] = size(obj.img.fullField);
            dff = zeros(T/2,N);
            for i=1:N
                temp = downsample(obj.img.fullField(:,i),2);                
                temp = zscore(temp ./ median(temp) - 1.0);
                dff(:,i) = temp';
            end
        end
            
        function [allDff, allDeconv, allSpikes] = Deconvolve(obj, traces, dt, useTime)
            % deconvolve and neuropil subtract
            warning ('off','all');
            if nargin < 3
                dt = 1/15.;
            end
            if nargin < 4
                useTime = 1:size(traces,2);
            end
            %traces = traces(:, useTime, :);
            allDff = zeros(size(traces));
            allDeconv = zeros(size(traces));
            allSpikes = zeros(size(traces));            
            [K, N, T] = size(traces);
            parfor_progress(K); % Initialize 
            parfor kk=1:K
                currTrace = squeeze(traces(kk,:,:));                
                [N, T] = size(currTrace);
                                                                
                flattenedTrace = reshape(currTrace, [size(currTrace,1)*size(currTrace,2), 1]);                

                
                % not sure what most of these really mean...               
                options = struct();
                options.dt = dt; 
                options.MaxIter = 10.;
                options.MaxInerIter = 50;
                options.TauStd = [0.2, 2];
                options.default_g = [0.6,0.9];
                options.method = 'cvx';
                b = [];
                c1 = [];
                g = [];
                sn = [];
                
                
                flattenedTrace(isnan(flattenedTrace)) = 0;
                flattenedTrace(isinf(flattenedTrace)) = 0;
                %
                %[c,~,~,~,~,sp] = constrained_foopsi(flattenedTrace,b,c1,g,sn,options);  
                try
                    [c,~,~,~,~,sp] = MCEM_foopsi(flattenedTrace,b,c1,g,sn,options);
                    %

                    allDeconv(kk,:,:) = reshape(c, [N T]);
                    allSpikes(kk,:,:) = reshape(sp, [N T]);      
                    allDff(kk,:,:) = reshape(flattenedTrace,[N T]);                
                catch
                    disp(['Could not deconvolve trial: ', num2str(kk)]);
                end
                parfor_progress;
            end
            parfor_progress(0);
        end
            
        function m = UnconcatenateTimeseries(obj, tseries, startPos, totalLen)
            % this could probably be done with reshape but whatever...            
            for i=1:size(tseries,3)
            end
        end
        
        function r = GetPerTrialLickRate(obj, binSize, maxTime)

            events = GetEventTimes(obj.behavior,'Port1In');

            r = EventRate(events, binSize, maxTime)';
        end
        
        function r = GetSmoothedPerTrialLickRate(obj, binSize, maxTime)
            sampleRate = 1/30.;
            events = GetEventTimes(obj.behavior,'Port1In');
            timePoints = sampleRate:sampleRate:maxTime;
            r = zeros(numel(timePoints), numel(events));
            for i=1:numel(events)
                currLicks = events{i};                
                lickITI = [0 1./diff(currLicks)];
                %for j=0:binSize:maxTime        
                for j=1:numel(timePoints)
                    t = timePoints(j);
                    currWin = [t t+binSize];
                    idx = find(currLicks >= currWin(1) & currLicks < currWin(2));
                    if numel(idx) > 0
                        r(j,i) = mean(lickITI(idx));
                    else
                        r(j,i) = 0; 
                    end%sum(currLicks >= currWin(1) & currLicks < currWin(2))/binSize;
                    %k = k + 1;
                end
            end

            %r = EventRate(events, binSize, maxTime)';
        end
        
        function aligned = AlignToEvent(obj, eventTimes)            
            
        end
        
        function [roc,trueRate,pval] = RocAnalyze(obj,timeWin)
            trialTypes = obj.behavior.TrialTypes;
            traces = obj.img.deconv(obj.img.taskLabels,timeWin,:);
            goTrials = trialTypes==1;
            nogoTrials = trialTypes==2;     
            roc = []; pval = [];
            for i=1:size(traces,1) % iterate over cells                
                currTraces = squeeze(traces(i,:,:));
                roc(i) = obj.doRoc(currTraces, trialTypes);
                % generate null distribution
                tempRoc = [];
                parfor kk=1:1000
                    tempTrials = trialTypes(randperm(numel(trialTypes)));
                    tempRoc(kk) = obj.doRoc(currTraces, tempTrials);
                end
                pval(i) = sum(roc(i) <= tempRoc)/numel(tempRoc);
            end
            trueRate = sum(trialTypes==1)/numel(trialTypes);
            trueRate = repmat(trueRate, 1, numel(roc));
        end
        
        function r = doRoc(obj, currTraces, trialTypes)
            dv = zeros(size(currTraces,2),1);
            for j=1:size(currTraces,2) % iterate over trials
                currOtherTrials = trialTypes(1:end ~= j);
                tempTraces = currTraces(:,1:end ~= j);
                meanGo = squeeze(mean(tempTraces(:,currOtherTrials == 1),2));
                meanNoGo = squeeze(mean(tempTraces(:,currOtherTrials==2),2));
                dv(j) = squeeze(currTraces(:,j))' * (meanGo-meanNoGo);
            end
            [~,~,~,r] = perfcurve(trialTypes, dv, 1);
        end
        
        function [correct, chance] = PopDecodeVariable(obj, labels, timeWin, classifierType,justActive)
            if nargin < 4
                classifierType = 'svm';
            end
            % data is Ntrial length vector
            traces = obj.img.deconv;
            if justActive
                traces = traces(obj.img.taskLabels,:,:);
            end
            if size(traces,1) > 1
                perTrialAvg = squeeze(mean(traces(:, timeWin, :),2));
                Nfold = 5; % fold for cross validation   
                performance = [];
                Ntrial = numel(labels);
                indices = crossvalind('Kfold',Ntrial,Nfold);
                chance = [];
                for i=1:Nfold
                    testIdx = (indices == i); trainIdx = ~testIdx;       
                    testLabels = labels(testIdx)';
                    trainLabels = labels(trainIdx)';
                    testData = perTrialAvg(:, testIdx)';
                    trainData = perTrialAvg(:, trainIdx)';

                    switch classifierType
                        case 'fisher'
                            wFisher = fisherLdaFit(trainData, trainLabels);
                            prediction = ((wFisher' * testData') > 0 + 1)';
                        case 'svm'
                            model = svmtrain(trainData, trainLabels);
                            prediction = svmclassify(model, testData);
                    end            
                    performance(end+1) = sum(prediction == testLabels)/numel(testLabels);
                    chance(end+1) = sum(testLabels==1)/numel(testLabels);
                end
                correct = mean(performance);
                chance = mean(chance);
            else
                correct = 0;
                chance = 0.5;
            end
            
        end
        
        function m = GetConcatenatedTimeseries(obj, tseries, startPos)
            if nargin < 3
                startPos = 1;
            end
            % assumes tseries in nCell x nTime x nTrial matrix
            m = [];
            for i=1:size(tseries,3)
                m = [m squeeze(tseries(:,startPos:end,i))];
            end
        end
        
        function r = MeanLicking(obj, trials, binSize, maxTime)
            r = EventRate(GetEventTimes(obj.behavior,'Port1In'), binSize, maxTime)';            
            if sum(trials)>1
                r = mean(r(trials,:));
            else
                r = r(trials,:);
            end
        end
        
        function c = LickingCorr(obj, trials, type, doMax, doZscore)
            % returns correlation of each cell with licking for trials
            frames = size(obj.img.deconv,2);
            maxTime = frames * obj.dt;
            
            r = obj.MeanLicking(trials, obj.dt, maxTime);
            if sum(r) == 0
                disp('here');
            end
            m = obj.MeanPerCellActivityForTrials(trials,type, doMax, doZscore);
            c = [];
            
            for i=1:size(m,1)
                c(i) = corr(m(i,:)', r(1:frames)');
            end 
            
        end
        
        function scaledCoords = GetCellLocs(obj)
            % get cell locations in global coordinate frame
            scalingFactor = 38/256;            
            if max(obj.img.cellCentroids(:)) < 256
                offset = -scalingFactor*128;
            else
                offset = -scalingFactor*256;
            end
            
             % hard coded for now = 6.7368 um/pixel in widefield
            goodCentroids = obj.img.cellCentroids(obj.goodTraces,:);             
            Ncell = size(goodCentroids,1);
            scaledCoords = zeros(Ncell,2);
            
            for i=1:Ncell
                scaledCoords(i,:) =  goodCentroids(i,:) * scalingFactor + offset;
            end            
        end
        
        function m = MeanPerCellActivityForTrials(obj,trials, type, doMax, doZscore)
            % returns mean activity for each trials per cell
            %m = obj.GetConcatenatedTimeseries(obj.img.traces, 30); % exclude first second
            %dff = obj.GetConcatenatedTimeseries(obj.img.neuropilTraces,30);
            if nargin < 3
                type = 'dff';
            end
            if nargin < 4
                doMax = false;
            end
            if strcmp(type, 'spikes')
                dff = obj.img.spikes;            
                for i=1:size(dff,1)
                    for j=1:size(dff,3)
                        dff(i,:,j) = smooth(dff(i,:,j),1/obj.dt,'moving');
                    end
                end
            elseif strcmp(type, 'deconv')                                
                dff = obj.img.deconv;
            else strcmp(type,'dff')
                dff = obj.img.traces;
            end
            if doZscore
                for i=1:size(dff,1)
                    for j=1:size(dff,3);
                        dff(i,:,j) = zscore(squeeze(dff(i,:,j)));
                    end
                end
            end
            
            m = squeeze(mean(dff(:,:,trials),3));
            if doMax
                for i=1:size(m,1)
                    m(i,:) = m(i,:)/max(m(i,:));
                end
            end
        end
        
        function m = MeanFullFieldActivityForTrials(obj,trials)
            m = mean(obj.img.fullField(:, trials),2);
        end
    end        
end

