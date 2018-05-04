function [zScore] = nanzscore(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if any(isnan(X(:)))
xmu=nanmean(X);
xsigma=nanstd(X);
zScore=(X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
else
[zScore,xmu,xsigma]=zscore(X);
end
end

