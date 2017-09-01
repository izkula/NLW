close all
clear all 

pathname='F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';

% calculate distances between the shock neurons or reward neurons, and
% randomly chosen neurons

% dates={'20160429','20160429','20160429','20160429','20160528'};
% filenames={'GRIN_pfcnac1_tonerewardshock_day5','GRIN_pfcnac2_tonerewardshock_day5','GRIN_pfcnac4_tonerewardshock_day5','GRIN_pfcnac5_tonerewardshock_day5','GRIN_pfcnac3_tonerewardshock_day4'};

dates={'20160529','20160529','20160730','20160730','20160731'};
filenames={'GRIN_pfcvta12_tonerewardshock_day5','GRIN_pfcvta15_tonerewardshock_day5','GRIN_pfcvta1_tonerewardshock_day4','GRIN_pfcvta2_tonerewardshock_day5','GRIN_pfcvta4_tonerewardshock_day6'};


dist_shock_all=[];
for z=1:length(dates)
    date=dates{z};
    filename=filenames{z};
    
    load(strcat(pathname,date,'/timecourses/',filename,'_cellmasks.mat'),'segcentroid');
    load(strcat(pathname,date,'/regression/',filename,'_shockonly.mat'));
    
    shock_centroids=segcentroid(shock_only,:);
    dist_shock=dist(shock_centroids');
    % distances in m/l direction
    %     dist_shock=dist(shock_centroids(:,2)');
    dist_shock=triu(dist_shock,1);
    dist_shock=dist_shock(dist_shock~=0);
    %dist_shock_all=[dist_shock_all; dist_shock];
    dist_shock_all=[dist_shock_all mean(dist_shock)];
    
    choose_shock(z)=length(shock_only);
    
    dist_nonmod_all=[];
    for p=1:1000
        nonmod_centroids=segcentroid;
        choose_nonmod=randperm(length(nonmod_centroids),choose_shock(z));
        nonmod_centroids=nonmod_centroids(choose_nonmod,:);
        
        dist_nonmod=dist(nonmod_centroids');
        % distances in m/l direction
        %         dist_nonmod=dist(nonmod_centroids(:,2)');
        dist_nonmod=triu(dist_nonmod,1);
        dist_nonmod=dist_nonmod(dist_nonmod~=0);
        dist_nonmod_all=[dist_nonmod_all; mean(dist_nonmod)];
        p
    end
    pval(z)=length(find(dist_nonmod_all<mean(dist_shock)))/1000
end