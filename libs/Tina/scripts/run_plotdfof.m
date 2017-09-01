function run_plotdfof(frametime,pathname,filedate,fn)

currpath=strcat(pathname,filedate);
load(strcat(currpath,'/timecourses/',fn,'_dfof.mat'),'dfof','dfof_sub','dfof_filt','dfof_deconv');

time=[0:frametime:(length(dfof)-1)*frametime];

figure(203);
for a=1:size(dfof,1)
     plot(time,dfof(a,:)+(a-1)*5);
%     plot(dfof(a,:)+(a-1)*5);
    hold on
end
xlabel('Time (s');
ylabel('dF/F');

% for a=1:size(dfof,1)
%     figure(200);
%     plot(time,dfof(a,:),'linewidth',2.5,'color',[25/255,25/255,112/255]);
%     hold on
%     plot(time,dfof_sub(a,:),'linewidth',2,'color',[0,0,205/255]);
%     hold on
%     plot(time,dfof_filt(a,:),'linewidth',1.5,'color',[65/255,105/255,225/255]);
%     hold on
%     %plot(time,dfof_deconv(a,:),'linewidth',1,'color',[135/255,206/255,235/255]);
%     
%     title(strcat('Cell ',num2str(a)));
%     hold off
%     
%     pause;
% end