function fish = inference(N, seed, knots_setting, est)
%Construct approximate confidence interval of r
%   approximate fisher information through E(grad^2)
if nargin < 4
    knots_setting = "quantile";
end

setting = 3;
load(strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat'), 'x', 'time', 'delta')

% desgin model
k = 4; % order 4 cubic spline
l = ceil(N^(1/7));% l polynomial pieces
[beta, ~, H] = coxphfit(x, time, 'Censoring', 1-delta);
[~, bin] = histc(time, H(:,1));
bin(bin == 0) = 1;
x_coef = exp(x * beta);
temp = x_coef .* H(bin,2);
if knots_setting == "quantile"
    knots = augknt([0, quantile(temp, (1:l-1)/l), max(temp)],k);
elseif knots_setting == "equal"
    knots = augknt(linspace(0, 2 * max(temp),l+1), k);
end
p = size(x,2);

forward = true;
beta = est(1:p);
theta = est(1+p:end);
[~, grad_beta] = objective_func_beta(beta, x, time, delta, theta, ...
    knots, k, true);
[~, grad_sieve] = objective_func_sieve(theta, x, time, delta, beta, ...
    knots, k, true, forward);
grad = [grad_beta grad_sieve];


grad(:, ~any(grad, 1)) = [];
% approximate fisher information
fish = inv((grad.' * grad));
end

