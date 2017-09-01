t_1=find(time_trial>=1,1,'first');
t_3=find(time_trial>=3,1,'first');

tina=mtrials_shock_press(shock_only,:);
ind=[];
for a=1:length(shock_only)
    %tina(a,:)=tina(a,:)./max(tina(a,:));
    [~,ind(a)]=max(tina(a,:));
    tina(a,:)=smooth(tina(a,:),4);
end

order=[ind' [1:length(shock_only)]'];
order=sortrows(order,1);

tina_sort=[];
for a=1:length(order)
    tina_sort(a,:)=tina(order(a,2),:);
end


figure;
imagesc(flipud(tina_sort));
colormap(jet);
colorbar
box off;
%caxis([0 5]);
hold on
plot([t_1 t_1],[0 length(shock_only)+1],'-w','linewidth',2)
plot([t_3 t_3],[0 length(shock_only)+1],'-w','linewidth',2)

%%
tina=mtrials_reward_press(shock_only,:);
for a=1:length(shock_only)
    %tina(a,:)=tina(a,:)./max(tina(a,:));
    tina(a,:)=smooth(tina(a,:),4);
end

tina_sort=[];
for a=1:length(order)
    tina_sort(a,:)=tina(order(a,2),:);
end

figure;
imagesc(flipud(tina_sort));
colormap(jet);
colorbar
box off;
%caxis([0 5]);
hold on
plot([t_1 t_1],[0 length(shock_only)+1],'-w','linewidth',2)
plot([t_3 t_3],[0 length(shock_only)+1],'-w','linewidth',2)

%%

tina=mtrials_reward_press(reward_only,:);
ind=[];
for a=1:length(reward_only)
    %tina(a,:)=tina(a,:)./max(tina(a,:));
    [~,ind(a)]=max(tina(a,:));
    tina(a,:)=smooth(tina(a,:),4);
end

order=[ind' [1:length(reward_only)]'];
order=sortrows(order,1);

tina_sort=[];
for a=1:length(order)
    tina_sort(a,:)=tina(order(a,2),:);
end


figure;
imagesc(flipud(tina_sort));
colormap(jet);
colorbar
box off;
%caxis([0 5]);
hold on
plot([t_1 t_1],[0 length(reward_only)+1],'-w','linewidth',2)
plot([t_3 t_3],[0 length(reward_only)+1],'-w','linewidth',2)

%%
tina=mtrials_shock_press(reward_only,:);
for a=1:length(reward_only)
    %tina(a,:)=tina(a,:)./max(tina(a,:));
    tina(a,:)=smooth(tina(a,:),4);
end

tina_sort=[];
for a=1:length(order)
    tina_sort(a,:)=tina(order(a,2),:);
end

figure;
imagesc(flipud(tina_sort));
colormap(jet);
colorbar
box off;
%caxis([0 5]);
hold on
plot([t_1 t_1],[0 length(reward_only)+1],'-w','linewidth',2)
plot([t_3 t_3],[0 length(reward_only)+1],'-w','linewidth',2)