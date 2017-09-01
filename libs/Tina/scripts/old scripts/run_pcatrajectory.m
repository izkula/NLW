function run_pcatrajectory(pathname,filedate,fn,frametime)

currpath=strcat(pathname,filedate);

% load lever, start_inds, dfof
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof_sub');
load(strcat(currpath,'/results/',fn,'.mat'));

dfof=dfof_sub;
time=[0:frametime:(length(dfof)-1)*frametime];

[F,pcs,eigs,var,pvar]=pca(dfof);

figure(500);
plot(pcs(:,1));
hold on
plot(pcs(:,2),'r');
plot(pcs(:,3),'g');
legend('PC1','PC2','PC3');
figure(501);
plot(pvar);





