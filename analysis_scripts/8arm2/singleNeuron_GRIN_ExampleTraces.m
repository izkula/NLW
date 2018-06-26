load('/home/svesuna/2presults/8arm/m864_8arm/neuron.mat')
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 2.5 2]; 

for i = 1:12
plot(C_raw(i,:)+i*50,'k');
hold on
plot(C(i,:)+ i*50,'r','LineWidth',.1);


end

%Make X Axis Correctly
xlabel('Time (s)')
x = (1:size(C,2))*(1/31);
x = round(x(1:31*120:end));
xticks = 1:size(C,2);
xticks = xticks(1:31*120:end);

set(gca,'XTick',xticks,'XTickLabel',x,'FontSize',5)

axis([0 1200*31 -inf inf])

saveas(h1,'~/Dropbox/2p_Claustrum_Shared/2p/Results/8arm/Figures/ExampleTraces', 'pdf')