function [ coeff, stats, opt_lim, log_X, y, mdl] = find_semipowerlawfit_with_grid_search(X, y, xlimits, grid_lim)

%modification to define all outputs made by Sam Lassers 9/1/21
coeff=[];
stats=[];
opt_lim=[];
log_X=[];
% y=[];

% Grid search xlim to find best r-squared
if ~exist('grid_lim','var')
    grid_lim = [0.5, 0.5]; %percentage value of changing x_start and x_stop
end
nbins = 50;

lim_grid = logspace(log10(xlimits(1)), log10(xlimits(2)), nbins);
% lim_grid = X;
x_start_vec = lim_grid(lim_grid < (1+grid_lim(1))*xlimits(1));
x_stop_vec = fliplr(lim_grid(lim_grid > (1-grid_lim(1))*xlimits(2)));

rsq_vec = []; start_vec = []; stop_vec = [];
for startn = x_start_vec
 for stopn = x_stop_vec
     [~, stats] = find_semipowerlawfit_using_linear_regression_coeff(X, y, [startn, stopn]);
     rsq_vec = [rsq_vec; stats.Rsquared];
     start_vec = [start_vec; startn];
     stop_vec = [stop_vec; stopn];
 end
end

opt_rowi = find(rsq_vec == max(rsq_vec),1);
opt_lim = [start_vec(opt_rowi), stop_vec(opt_rowi)];
%modification by Sam Lassers 9/1/21 in case of low population vector
%causing optimal limit to not be found
if ~isempty(opt_lim)
    [ coeff, stats, log_X, y, mdl] = find_semipowerlawfit_using_linear_regression_coeff( X, y, opt_lim );
else
    warning('Optimal limits for regression not found. Regression will not be plotted.')
end
end
 
         