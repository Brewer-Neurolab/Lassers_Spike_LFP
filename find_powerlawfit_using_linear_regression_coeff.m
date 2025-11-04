function [ coeff, stats, log_X, log_y ] = find_powerlawfit_using_linear_regression_coeff( X, y, xlim)
%Find linear regression coefficients between X and y in range [xlim(1),
%xlim(2)]
if ~exist('xlim','var')
    xlim = [-inf, inf];
end
% removing zeros
zI = X <= 0 | y <=0;
X(zI) =[]; y(zI) = [];
% applying limits specified in xlim
nzI =  X > xlim(1) & X < xlim(2);
log_X = log10(X(nzI)); log_y = log10(y(nzI));
mdl = fitlm(log_X, log_y);
coeff = mdl.Coefficients.Estimate;
stats = struct;
stats.Rsquared = mdl.Rsquared.Ordinary;
stats.pValue = mdl.Coefficients.pValue;
stats.coeffNames = {'Intercept','Slope'};

% stats.
end

