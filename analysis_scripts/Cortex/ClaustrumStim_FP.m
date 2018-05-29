clear; clc; 
init_samz
nidaqSampleRate = 5000;
downSampleFactor = 500

f = {'/media/svesuna/Sam_2p_Proc_Data_B/20180417/m845_stim_fp/nidaq.mat'
    '/media/svesuna/Sam_2p_Proc_Data_B/20180417/m846_stim_fp/nidaq.mat'
    '/media/svesuna/Sam_2p_Proc_Data_B/20180417/m861_stim_fp/nidaq.mat'
    }


for ii = 1:numel(f)
    [t,d] = LoadNidaqOutput(f{ii},5);
    
    
     %%%% EXTRACT RUNNING DATA %%%%%
    runTicks = d(5,:) < -0.1;
    tickBin = nidaqSampleRate / 10; %for a bin of 100 ms
    runCounts = sum(reshape(runTicks, tickBin, []));
    %runCountsSmooth = smooth(runCounts, 50); 
   
    
    %%%%%% EXTRACT FP SIGNAL %%%%%
    fp = smooth(d(3,:) - d(4,:),1000);
    fp = mean(reshape(fp, tickBin,[]));
    
    
    tStarts = find(diff(d(2,:) > 3.99999999));
    tStarts = tStarts(1:2:end) / downSampleFactor;

    nSamplesBefore = 4* nidaqSampleRate / downSampleFactor
    nSamplesStim = 2* nidaqSampleRate / downSampleFactor;
    nSamplesAfter = 4* nidaqSampleRate / downSampleFactor;
    nSamplesTotal = (nSamplesBefore + nSamplesAfter + nSamplesStim);
    

    fp_t = zeros(10,nSamplesTotal); fp_stats = []; run_stats = [];
    for i = 1:10
        fp_t(i,:) = fp(tStarts(i) - nSamplesBefore +1 : tStarts(i) + nSamplesStim + nSamplesAfter); 
        
        
        for j = 1:5
            x1= (j-1)*20+1;
            x2= j*20;
             fp_stats(i,j) = mean(fp_t(i,x1:x2)) ;
        end
        
        %run_t = runCounts(tStarts(i) - nSamplesBefore : tStarts(i) + nSamplesStim + nSamplesAfter-1);
       % run_stats(i,1:5) = [mean(fp_t(i,1:10)), mean(fp_t(i,10:20)), mean(fp_t(i,20:30)), mean(fp_t(i,30:40)), mean(fp_t(i,40:50))] ;

        
    end

    

    
    
%%%% CALCULATE SOME SHIT %%%%
%calculate mean
fp_t = fp_t - mean(fp_t(:,1:nSamplesBefore),2);

%calculate fp statistics
[P, ANOVATAB,STATS] = anova1(fp_stats(:,2:4));
C = multcompare(STATS);

%calculate running correlation
[rho, p] = corr(fp', runCounts')


%Plot
h1 = figure; h1.Units = 'inches'; h1.Position = [3 3 3 3];

%plot 1
subplot(2,2,[1 3])
%plot(fp_t','k','LineWidth',.1)
rectangle('Position',[nSamplesBefore min(mean(fp_t)) nSamplesStim max(mean(fp_t) - min(mean(fp_t))) ],'FaceColor','r')

hold on; plot(mean(fp_t),'k','LineWidth',3)

secTick = 1; %seconds for ticks
x = (1:nSamplesTotal)*1/10;
x = x(1:10*secTick:end);
xticks = 1:nSamplesTotal;
xticks = xticks(1:10*secTick:end);

set(gca,'XTick',xticks,'XTickLabel',round(x), 'FontSize', 6)


xlabel 'Time (s)'
ylabel 'dF (F - baseline)'

%plot 2
hold off
subplot(2,2,[2])

errorb(mean(fp_stats(:,2:4)),sem(fp_stats(:,2:4)))
hold on; bar(mean(fp_stats(:,2:4)));
plot(fp_stats(:,2:4)','k','LineWidth',.1)
set(gca,'XTick',[1:3],'XTickLabel',{'Off' 'On' 'Off'}, 'FontSize', 6)

s = split(f{ii},'/');
subject = s{6};

%title
title([' p = ' num2str(C(1,end))])

%plot 3
hold on
subplot(2,2,4); 

x1 = fp - min(fp);
x2 = x1 ./ max(x1);
y1 = runCounts - min(runCounts);
y2 = y1 ./ max(y1);


plot(x2,'k');  
plot(y2,'r');
title(['rho =' num2str(rho)])

%save out
saveas(h1,fullfile('~/Dropbox/2p_Claustrum_Shared/2p/Results/Cortex/Figures/',subject), 'pdf')




end
