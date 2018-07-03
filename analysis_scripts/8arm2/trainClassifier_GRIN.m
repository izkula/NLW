clear;
load('/home/svesuna/2presults/8arm2/allNeuronsTrial.mat')

%%

%PKS = cell(numel(CT_S),1)
%PKS_Spikes = cell(numel(CT_S),1)

firstFrame = 124 - 3*31
lastFrame = nFramesTrial - 7*31
counter = 1;
for i = 1:numel(CT_S)
    for j = 1:size(CT_S{i},1)
        %extract calcium peak times
        for k = 1:size(CT_S{i},3)
            %PKS{i}(j,k) = length(findpeaks(CT_S{i}(j,firstFrame:lastFrame,k)));
            
           % PKS_Spikes{i}(j,k) = sum(ST(counter,firstFrame:lastFrame,k));
            Ca_avg{i}(j,k) = mean(CT_S{i}(j,firstFrame:lastFrame,k),2)
        end
        counter = counter + 1
    end
    i
end

%%

%%
%INDIVIDUAL NEURON PLOTTING FOR DIFFERENT TRIAL TYPES
trialTypes = [];
trialTypes = repmat([1 4 2 0 5 3 6 0],1,8)

 acc_base = [];

 j = 5
 temp = [];
    temp = [Ca_avg{j}; trialTypes];

    data = []; counter = 1
    for i = 1:size(temp,2)
        if temp(end,i) ~=0
            data(:,counter) = temp(:,i);
            counter = counter+1 
        end
    end

%%

acc_base = []; acc_shuffle = [];



parfor ii = 1:10
[t, acc_base(ii)] = trainClassifier(data);
end
mean(acc_base)

D = {}
for i = 1:1000
    data(end,:) = data(end,randperm(size(data,2)));
   
    D{i} = data;

end

parfor i = 1:1000
   [t, acc_shuffle(i)] = trainClassifier(D{i});
end
mean(acc_shuffle)

accuracySVM{j} = [acc_base, acc_shuffle]




%%
h1 = figure; h1.Units = 'inches'; h1.Position = [1 1 2 3];
axis([ 0.5 5.5 0 1])
counter = 1
for i = [1 4 5 6 7]
    try
    hold on
    %errorb(counter, mean(accuracySVM{i}(1:10)), std(accuracySVM{i}(1:10)'))
    plot(counter, mean(accuracySVM{i}(1:10)),'ob', 'MarkerSize', 3, 'MarkerFaceColor', 'g','MarkerEdgeColor','g')
    errorb(counter, mean(accuracySVM{i}(11:end)), std(accuracySVM{i}(11:end)'),'LineWidth',1)
        counter = counter+1
    end
end


set(gca,'FontSize',5)
ylabel('Classification Accuray')
xlabel('Mouse ID')
saveas(h1, [fig_dir '/Classifier'],'pdf')

 
%% signifiance


for j = [1 4 7]
    x = mean(accuracySVM{j}(1:10))
    
    y = sort(accuracySVM{j}(11:end));
    
    z = sum(y>x) / 1000
    
    
    
end

