load('~/2presults/8arm2/allNeuronsTrial.mat')
%%
rho = {}; p = {}; rho_proc = {};
for i = 1:numel(SSS)
    
    s = SSS{i}';
    
    s(isnan(s)) = 0;
    
    [rho{i}, p{i}] = corr(s,'rows','all');
    
    
    rho_proc{i} = reshape(triu(rho{i},1),1,[]);
    
    rho_proc{i}(rho_proc{i} == 0) = NaN;
   
end
