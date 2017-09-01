classdef MultiSession2P
    %DATASET2P collection of 2P collected from the same mouse, from the
    %same day
    %Includes optional calibration image and coordinates for each dataset
    
    properties
        datasets = {};
        roiCoords = [];
        registeredCoords = [];
        dates = {};
        behavior = [];
        calPath = '~/Dropbox/2p_cal';
        Ncells = 0;
        Nt = 0;
        dt = 0;
        
        oegDataPath = '';
        oegBehaviorPath = '';
        oegCoords = [];
        oegTraces = [];
        oegRef = [];
        
        refPointFront = [];
        refPointBack = [];
        maxT = 0;
        
        calImage = containers.Map();
        refImage = [];
         
        %coregPixelLocs = [];
    end
    
    methods
        function obj = MultiSession2P()
            obj.calImage = containers.Map();
        end

        function img = LoadCalImage(obj, mouseName, date)
            disp(fullfile(obj.calPath, date, mouseName, [mouseName '_MMStack.ome.tif']));
           img = imread(fullfile(obj.calPath, date, mouseName, [mouseName '_MMStack.ome.tif']));
        end
        
        
            
        function obj = LoadData(obj, mouseName, protocolName, date, session, loc, dt,maxT)                       
            obj.dt = dt;
            obj.maxT = maxT;
            obj.dates{end+1} = date;                        
            obj.calImage(date) = obj.LoadCalImage(mouseName, date);                        
            obj.roiCoords = [obj.roiCoords; loc];
            obj.datasets{end+1} = Dataset2P(mouseName, protocolName, date, session, dt,maxT);
        end
        
        function obj = LoadDataFromMAT(obj, mouseName, protocolName, date, session, loc, dt, maxT)
            obj.dt = dt;
            obj.maxT = maxT;
            obj.dates{end+1} = date;
            obj.calImage(date) = obj.LoadCalImage(mouseName, date);
            obj.roiCoords = [obj.roiCoords; loc];
            temp = load(fullfile('~/Dropbox/oeg_results/2p', mouseName, protocolName, date, ['Session' num2str(session)], 'deconvolved.mat'));
            obj.datasets{end+1} = temp.dataset;
        end
        
        function obj = SetRefPoints(obj)
            figure, imagesc(obj.refImage);
            obj.refPointFront = ginput(1);
            obj.refPointBack = ginput(1);
        end
        
        function dists = GetROIDists(obj)
            pts = obj.registeredCoords;            
            refPts = [obj.refPointFront; pts];
            dists = squareform(pdist(refPts,'euclidean'));
            dists = dists(1,2:end);
        end
        
        function DrawROILocs(obj, N)
            if nargin  < 2
                N = size(obj.registeredCoords,1);
            end
            if N > size(obj.registeredCoords,1)
                N = size(obj.registeredCoords,1);
            end
            locs = obj.registeredCoords(1:N,:);
            imagesc(obj.GetCalImage());  hold on; 
            for i=1:size(locs,1)
                XX = locs;%GetPixelLoc(locs(i,1), locs(i,2));
                %plot(XX(:,1), XX(:,2), 'om'); 
                rectangle('position', [XX(i,1)-19, XX(i,2)-19, 38, 38], 'EdgeColor',[1 0 0]); hold on;
                text(XX(i,1), XX(i,2), num2str(i), 'color', 'r');
                %x = round(locs(i,1)/13)+x0; y = round(locs(i,2)/13)+y0;
                %temp((x-19):(x+19), (y-19):(y+19)) = 0;
                hold on;
                
            end            
            colormap(gray(256));
        end
        
        function ColorROILocs(obj, vals)
            % vals is of length obj.registeredCoords
            cols = bluewhitered(200);
            field = obj.registeredCoords;
            %imagesc(obj.GetCalImage());  colormap(gray); hold on; 
            for i=1:numel(vals)
                if vals(i) > 0
                    %plot(field(i,1), field(i,2), '.','markersize',40,'color',cols(round(vals(i)),:)); hold on;
                    scatter(field(:,1),field(:,2),40, vals,'filled'); colormap(bluewhitered); caxis([-1 1]);
                end
            end

        end
        % stuff relating to matching with OEG dataset        
        function obj = SetOEGPath(obj, mouseName, protocol, date, session)
            obj.oegDataPath = fullfile('~/Dropbox/oeg_results/rcgc/', [mouseName '_' protocol '_' date '_Session' num2str(session)], 'avgVids_zscore.mat');
            obj.oegBehaviorPath = fullfile('~/Dropbox/oeg_results/rcgc/', [mouseName '_' protocol '_' date '_Session' num2str(session)], 'bpod.mat');
        end
                
        function [roc, trueRates,pvals] = GetAllROC(obj,timeWin)
            roc = []; trueRates = []; pvals = [];
            for currIdx=1:numel(obj.datasets)
                [r,t,p] = obj.datasets{currIdx}.RocAnalyze(timeWin);
                roc = [roc r];
                trueRates = [trueRates t];
                pvals = [pvals p];
            end
        end
        
        function [go, nogo] = GetAllLickingCorrs(obj,doZscore)
            % get per cell licking correlation
            go = []; nogo = [];
            for currIdx=1:numel(obj.datasets)
                currLabels = obj.datasets{currIdx}.img.taskLabels;    
                goTrials = obj.datasets{currIdx}.behavior.TrialTypes==1 & obj.datasets{currIdx}.success == 1;
                nogoTrials = obj.datasets{currIdx}.behavior.TrialTypes==2 & obj.datasets{currIdx}.success == 1;
                r = obj.datasets{currIdx}.GetPerTrialLickRate(0.033, 8);
                
                if sum(goTrials) > 1
                    rGo = mean(r(goTrials,:))';
                else
                    rGo = r(goTrials,:)';
                end
                if sum(nogoTrials) > 1
                    rNoGo = mean(r(nogoTrials,:))';
                else
                    rNoGo = r(nogoTrials,:)';
                end
                
                currCellsGO = obj.datasets{currIdx}.MeanPerCellActivityForTrials(goTrials,'deconv', false, doZscore);    
                currCellsNOGO = obj.datasets{currIdx}.MeanPerCellActivityForTrials(nogoTrials,'deconv', false,doZscore);    
                %currCellsGO = currCellsGO(currLabels,:);
                %currCellsNOGO = currCellsNOGO(currLabels,:);
                   
                cGo = [];
                cNoGo = [];
                
                    for i=1:size(currCellsGO)
                        cGo(i) = corr(rGo(1:obj.maxT), currCellsGO(i,:)');
                    end
                
                    for i=1:size(currCellsNOGO)
                        cNoGo(i) = corr(rNoGo(1:obj.maxT), currCellsNOGO(i,:)');
                    end
                
                go = [go cGo];
                nogo = [nogo cNoGo];
                
                %if size(lickGo,1) > 0
                %    cellSource = [cellSource; repmat(k, size(lickGo,1), 1)];
                %end
            end         
            
        end
        
        function [go, nogo] = GetAllAvgLicking(obj)
            k = 1; go = []; nogo = [];
            binSize = obj.dt;
            maxTime = obj.maxT * binSize;
            for currIdx=1:numel(obj.datasets)
                lickGo = obj.datasets{currIdx}.MeanLicking(obj.datasets{currIdx}.behavior.TrialTypes==1 & obj.datasets{currIdx}.success == 1, binSize, maxTime);    
                lickNogo = obj.datasets{currIdx}.MeanLicking(obj.datasets{currIdx}.behavior.TrialTypes==2 & obj.datasets{currIdx}.success == 1, binSize, maxTime);    
                go = [go lickGo'];
                nogo = [nogo lickNogo'];
                %if size(lickGo,1) > 0
                %    cellSource = [cellSource; repmat(k, size(lickGo,1), 1)];
                    k = k + 1;
                %end
            end         
            go = go(1:obj.maxT,:); nogo = nogo(1:obj.maxT,:);
        end
               
        function good = GetAllGoodCells(obj)
            good = [];
            for currIdx=1:numel(obj.datasets)
                good = [good 1:obj.datasets{currIdx}.img.Ncells];
            end
        end
        
        function active = GetAllActiveCells(obj)
            active = [];
            for currIdx=1:numel(obj.datasets);
                active = [active; obj.datasets{currIdx}.img.taskLabels];
            end
        end
        
        function cellLocs = GetAllCellLocs(obj)
            % get the locations of  cells
            cellLocs = [];
            for currIdx=1:numel(obj.datasets);
                locs = obj.datasets{currIdx}.GetCellLocs();
                currCellLocs = locs + repmat(obj.registeredCoords(currIdx,:), [size(locs,1),1]);
                cellLocs = [cellLocs; currCellLocs];
            end
        end
        
        function p = GetAllPC(obj,nPC)            
            p = cell(numel(obj.datasets),1);
            for currIdx=1:numel(obj.datasets)
                [go,~] = obj.datasets{currIdx}.GetTrialTypeAvg(true,true);
                currPC = mypca(go',nPC);
                p{currIdx} = currPC;
            end
        end
        
        function [go, nogo, cellSource] = GetAllAvgCells(obj, doZscore, justActive)
            % gets average activity on hit and cr trials
            go = []; nogo = []; cellSource = [];
            %{
            for i=1:numel(obj.datasets)
                currData = obj.datasets{i};
                [currGo, currNogo] = currData.GetTrialTypeAvg(doZscore, justActive);
                cellSource = [cellSource; repmat(i, size(currGo,1), 1)];
                go = [go; currGo];
                nogo = [nogo; currNogo];
            end
            %}
            k = 1;
            for currIdx=1:numel(obj.datasets)
                currLabels = obj.datasets{currIdx}.img.taskLabels;    
                currCellsGO = obj.datasets{currIdx}.MeanPerCellActivityForTrials(obj.datasets{currIdx}.behavior.TrialTypes==1 & obj.datasets{currIdx}.success == 1,'deconv', false, doZscore);    
                currCellsNOGO = obj.datasets{currIdx}.MeanPerCellActivityForTrials(obj.datasets{currIdx}.behavior.TrialTypes==2 & obj.datasets{currIdx}.success == 1,'deconv', false,doZscore);    
                %currLabels = ones(size(currCellsGO,1),1);
                %currCellsGO = currCellsGO(currLabels,:);
                %currCellsNOGO = currCellsNOGO(currLabels,:);
                go = [go currCellsGO'];
                nogo = [nogo currCellsNOGO'];
                if size(currCellsGO,1) > 0
                    %cellSource = [cellSource; repmat(k, size(currCellsGO,1), 1)];
                    cellSource = [cellSource; repmat(currIdx, size(currCellsGO,1), 1)];
                    k = k + 1;
                end

            end

        end
        
        function [oegData, oegBehavior] = LoadOEGDataset(obj)
            disp('loading oeg data');
            oegData = load(obj.oegDataPath);
            oegBehavior = load(obj.oegBehaviorPath);            
        end
        
        function obj = ExtractOEGTraces(obj)
            [oegData, oegBehavior] = obj.LoadOEGDataset();
            go = []; goV = []; goN = oegData.nGo;
            nogo = []; nogoV = []; nogoN = oegData.nNogo;
            coords = obj.oegCoords;
            for i=1:size(coords,1)                
                go = [go squeeze(oegData.G_GO_z(coords(i,1),coords(i,2),:))];
                goV = [goV squeeze(oegData.G_GO_z_var(coords(i,1),coords(i,2),:))];
                nogo = [nogo squeeze(oegData.G_NOGO_z(coords(i,1),coords(i,2),:))];
                nogoV = [nogoV squeeze(oegData.G_NOGO_z_var(coords(i,1),coords(i,2),:))];
            end                   
            obj.oegTraces = struct();
            obj.oegTraces.go = go;
            obj.oegTraces.goV = goV;
            obj.oegTraces.nogo = nogo;
            obj.oegTraces.nogoV = nogoV;
        end
        
        function [c,idx] = GetSortedRegisteredCoords(obj)
            c = obj.registeredCoords;
            [~,idx] = sort(sqrt(c(:,2).^2 + c(:,1).^2));
            c = c(idx,:);
        end
        
        function im = GetCalImage(obj)
            k = keys(obj.calImage);
            im = obj.calImage(k{1});
        end
        
        function obj = CoregisterMultipleDays(obj,N)
            if nargin < 2
                N = 4;
            end
            pixelLocs = zeros(size(obj.roiCoords));
            for i=1:size(obj.roiCoords,1)
                pixelLocs(i,:) = GetPixelLoc(obj.roiCoords(i,1), obj.roiCoords(i,2));
            end
            days = keys(obj.calImage);
            obj.refImage = obj.calImage(days{1});
            
            disp('click on N points in each image');            
            figure(1), imagesc(obj.calImage(days{1})); colormap(gray(256));
            coordsRef = ginput(N);            
            regs = containers.Map();
            otherCoords = {};
            currDates = {};
            
            for i=2:numel(days)
                figure(1), imagesc(obj.calImage(days{i})); colormap(gray(256));
                otherCoords{end+1} = ginput(N);
                currDates{end+1} = days{i};
            end
            for i=1:numel(otherCoords)
                temp = otherCoords{i};
                regs(currDates{i}) = estimateGeometricTransform(temp, coordsRef, 'similarity');                
            end
            obj.registeredCoords = [];
            for i=1:size(obj.roiCoords,1)
                if strcmp(obj.dates{i}, days(1))
                    obj.registeredCoords = [obj.registeredCoords; pixelLocs(i,:)];
                else                   
                    currReg = regs(obj.dates{i});
                    temp = pixelLocs(i,:);
                    [x,y] = transformPointsForward(currReg, temp(1), temp(2));
                    obj.registeredCoords = [obj.registeredCoords; [x y]];
                end
            end
            %obj.coregPixelLocs = pixelLocs;  
        end
        
        function  obj = FindOEG2PCorrespondance(obj,N)
            % returns points in OEG vid that correspond to 2P locations
            % requires having loaded cal image and 2P dataset
            if nargin < 2
                N = 2;
            end
            [oegData,~] = obj.LoadOEGDataset();
            obj.oegRef = oegData.G_GO_z(:,:,1);
            disp('click on 2 points in each image to find correspondence');
            figure('position', [100 100 1024 512]);
            subplot(1,2,1);
            imagesc(obj.oegRef); colormap(gray(256)); axis square; hold on;
            gOEG = ginput(N);
            plot(gOEG(:,1), gOEG(:,2), 'r*');
            subplot(1,2,2);
            imagesc(obj.refImage); colormap(gray(256)); axis square; hold on;
            gCal = ginput(N);
            plot(gCal(:,1), gCal(:,2), 'r*');

            tform = estimateGeometricTransform(gCal, gOEG, 'similarity');            
            %inverseTform = estimateGeometricTransform(gOEG, gCal, 'similarity');            
            c = obj.GetSortedRegisteredCoords();
            [x,y] = transformPointsForward(tform, c(:,1), c(:,2));
            warpedCal = imwarp(obj.refImage, tform, 'OutputView', imref2d(size(obj.oegRef)));
            %[xR,yR] = transformPointsInverse(tform, pixelLocs(:,1), pixelLocs(:,2));
            figure(1); clf;
            subplot(1,2,1);    
            
            imagesc(obj.oegRef); colormap(gray(256)); hold on;% plot(x,y,'r*'); axis off;
            for i=1:numel(x)
                text(x(i), y(i), num2str(i), 'color','r');
            end
            axis off;
            subplot(1,2,2);            
            imagesc(warpedCal); colormap(gray(256)); hold on; %plot(pixelLocs(:,1), pixelLocs(:,2),'r*'); axis off;
            for i=1:numel(x)
                text(x(i), y(i), num2str(i), 'color','r');
            end            
            axis off;
            obj.oegCoords = round([x y]);
        end
         
        function [cells, session] = GetAllCells(obj)
            session = [];
            cells = [];
        end
    end
    
end

