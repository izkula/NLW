close all
clear all
pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';


dates={'20160429','20160429','20160429','20160429','20160528'};
filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

pressed=[];
missed=[];
nac_mice=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/regression/',filename,'_peaks.mat'));
%     peaks_missed(peaks_missed<0)=0;
%     peaks_pressed(peaks_pressed<0)=0;
    
    missed=[missed peaks_missed];
    pressed=[pressed peaks_pressed];
    temp=(missed-pressed)./(missed+pressed);
    nac_mice(z)=mean(temp(~isnan(temp)));
end

si_nac=(missed-pressed)./(missed+pressed);
%si_nac=(missed-pressed)./(pressed);
si_nac(isnan(si_nac))=[];

nac_missed=missed;
nac_pressed=pressed;

figure;
bar([1:2],[mean(nac_missed) mean(nac_pressed)])
hold on
errorbar([1:2],[mean(nac_missed) mean(nac_pressed)],[std(nac_missed)/sqrt(length(nac_missed)) std(nac_pressed)/sqrt(length(nac_pressed))])
signrank(nac_missed,nac_pressed)





dates={'20160529','20160529','20160730','20160730','20160731'};
filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};

pressed=[];
missed=[];
vta_mice=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};

    load(strcat(pathname,date,'/regression/',filename,'_peaks.mat'));
%     peaks_missed(peaks_missed<0)=0;
%     peaks_pressed(peaks_pressed<0)=0;
    
    missed=[missed peaks_missed];
    pressed=[pressed peaks_pressed];
    temp=(missed-pressed)./(missed+pressed);
    vta_mice(z)=mean(temp(~isnan(temp)));
end

si_vta=(missed-pressed)./(missed+pressed);
%si_vta=(missed-pressed)./(pressed);
si_vta(isnan(si_vta))=[];

vta_missed=missed;
vta_pressed=pressed;

figure;
bar([1:2],[mean(vta_missed) mean(vta_pressed)])
hold on
errorbar([1:2],[mean(vta_missed) mean(vta_pressed)],[std(vta_missed)/sqrt(length(vta_missed)) std(vta_pressed)/sqrt(length(vta_pressed))])
signrank(vta_missed,vta_pressed)

figure;
bar([1,2],[mean(si_nac),mean(si_vta)]);
hold on
errorbar([1,2],[mean(si_nac),mean(si_vta)],[std(si_nac)/sqrt(length(si_nac)),std(si_vta)/sqrt(length(si_vta))]);

p_nac=signtest(si_nac)
p_vta=signtest(si_vta)

% figure;
% bar([1,2],[mean(nac_mice),mean(vta_mice)]);
% hold on
% errorbar([1,2],[mean(nac_mice),mean(vta_mice)],[std(nac_mice)/sqrt(length(nac_mice)),std(vta_mice)/sqrt(length(vta_mice))]);