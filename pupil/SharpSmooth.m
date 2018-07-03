function smoothed = SharpSmooth(covariate, thresh, smoothFactor)
%%% Smooths the provided covariate, but sharply clamps the start and beginning of
%%% epochs, where epochs are nontrivial stretches of activity above
%%% the provided thresh. 

epochs = FindRunEpochs(covariate, thresh);
causalSmoothed = CausalSmooth(covariate, smoothFactor)'; %
smoothed = epochs.*causalSmoothed;
%smoothed = epochs;
