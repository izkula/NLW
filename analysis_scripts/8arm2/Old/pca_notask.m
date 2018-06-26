clear; clc

load('/home/svesuna/2presults/Running/m76_run/z0_neuron.mat')






























if isempty(find(isnan(C)))
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(C);
else
   C1 = fillmissing(C,'previous',2);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(C1);

end

figure
plot(cumsum(EXPLAINED))
sFactor = 50

%
figure
c = linspecer(length(COEFF))
h = colormapline(smooth(COEFF(:,1),sFactor), smooth(COEFF(:,2),sFactor), [], c)
set(h,'linewidth',3)

%%
plot3(smooth(COEFF(:,1),sFactor), smooth(COEFF(:,2),sFactor), smooth(COEFF(:,3),sFactor))