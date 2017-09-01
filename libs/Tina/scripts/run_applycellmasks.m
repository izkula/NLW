function run_applycellmasks(pathname,filedate,fn,frametime,deconvtauoff,slidingwind,percentfilter,lowpassfilter)

currpath=strcat(pathname,filedate);

outputdir=strcat(currpath,'/timecourses/');

% Get cell time series
load(strcat(currpath,'/timecourses/',fn,'_cellmasks.mat'));

subtractmean=0;
flims=[];
% cell_sig = CellsortApplyFilter_removenoise(strcat(currpath,fn,'_mc.tif'), cellmask, flims, subtractmean);
cell_sig = CellsortApplyFilter_removenoise(strcat(currpath,fn,'.tif'), cellmask, flims, subtractmean);


% Set parameters
[dfof,dfof_sub,dfof_filt,dfof_deconv] = GetTimeSeries(cell_sig,frametime,deconvtauoff,slidingwind,percentfilter,lowpassfilter);

save(strcat(outputdir,fn,'_dfof.mat'),'cell_sig','dfof','dfof_sub','dfof_filt','dfof_deconv');

function [dfof,dfof_sub,dfof_filt,dfof_deconv] = GetTimeSeries(cell_sig,frametime,deconvtauoff,slidingwind,percentfilter,lowpassfilter)

% Get dfof of cells
dfof=zeros(size(cell_sig,1),size(cell_sig,2));
% for a=1:size(cell_sig,1)
%     cell=cell_sig(a,:);
%     baseline=median(cell);
%     dfof(a,:)=(cell-baseline)/baseline;
% end
for a=1:size(cell_sig,1)
    cell=cell_sig(a,:);
    %baseline=median(cell);
    dfof(a,:)=zscore(cell);
end

% subtract nth percentile value over sliding window
wind=2.*round(slidingwind/frametime/2);
dfof_sub=zeros(size(dfof,1),size(dfof,2));
for a=1:size(dfof,1)
    subval=zeros(1,size(dfof,2));
    for b=1:size(dfof,2)
        if b<=wind/2
            subval(b)=prctile(dfof(a,1:b+wind/2),percentfilter);
        elseif b>wind/2 && b<=size(dfof,2)-wind/2
            subval(b)=prctile(dfof(a,b-wind/2:b+wind/2),percentfilter);
        elseif b>size(dfof,2)-wind/2
            subval(b)=prctile(dfof(a,b-wind/2:end),percentfilter);
        end
    end
    dfof_sub(a,:)=dfof(a,:)-subval;
    offset=median(dfof_sub(a,:));
    dfof_sub(a,:)=dfof_sub(a,:)-offset;
end

% low pass filter
fNorm = lowpassfilter / (1/frametime/2);
[b, a] = butter(2, fNorm, 'low');
dfof_filt=zeros(size(dfof,1),size(dfof,2));
for c=1:size(dfof,1)
    dfof_filt(c,:) = filtfilt(b, a, dfof_sub(c,:));
end

% Calculate deconvolved timeseries
dfof_deconv=zeros(size(dfof,1),size(dfof,2));
tdc=[0:0.01:0.5];
u=exp(-tdc./deconvtauoff);
u=u./sum(u);
dfof_filt_pad=[dfof_filt zeros(size(dfof_filt,1),length(u)-1)];
for d=1:size(dfof,1)
    dfof_deconv(d,:)=deconv(dfof_filt_pad(d,:),u);
end

display('Completed running applycellmasks.m');
