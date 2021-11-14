function fish = inference(N, seed, knots_setting, r)
%Construct approximate confidence interval of r
%   approximate fisher information through E(grad^2)

% load data
setting = 2;
if nargin < 4
    knots_setting = "quantile";
end
rho1 = 0; r1 = 1;
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
load(data_file, 'x', 'time', 'delta');

% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces scales with the sample size
temp = time(delta > 0);
if knots_setting == "quantile"
    knots = augknt([0, quantile(temp, (1:l-1)/l), max(time)],k);
elseif knots_setting == "equal"
    knots = augknt(linspace(0,max(time),l+1), k);
end
q = size(knots, 2) - k;

[~, grad] = objective_func(r, x, time, delta, q, knots, k, rho1, r1, true);

% grad(:, ~any(grad, 1)) = [];

% approximate fisher information
fish = inv((grad.' * grad));
end