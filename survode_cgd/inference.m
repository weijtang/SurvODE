function fish = inference(N, seed, setting, knots_setting, est)
%Construct approximate confidence interval of est
%   approximate fisher information through E(grad^2)

data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
load(data_file, 'x', 'time', 'delta');

[N, p] = size(x);
k0 = 4; % order 4 cubic spline
kq = 4; 
l0 = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces for alpha()
lq = ceil(N^(1/7)); % l polynomial pieces for q()

temp1 = time(delta > 0);
[beta, ~, H] = coxphfit(x, time, 'Censoring', 1-delta);
[~, bin] = histc(time, H(:,1));
bin(bin == 0) = 1;
x_coef = exp(x * beta);
temp2 = x_coef .* H(bin,2);

if knots_setting == "K1"
    knots_0 = augknt(linspace(0,max(time),l0+1), k0);
    knots_q = augknt(linspace(0, 2 * max(temp2),lq+1), kq);
elseif knots_setting == "K2"
    knots_0 = augknt(linspace(0,max(time),l0+1), k0);
    knots_q = augknt([0, quantile(temp2, (1:lq-1)/lq), max(temp2)],kq);
elseif knots_setting == "K3"
    knots_0 = augknt([0, quantile(temp1, (1:l0-1)/l0), max(time)],k0);
    knots_q = augknt(linspace(0, 2 * max(temp2),lq+1), kq);
elseif knots_setting == "K4"
    knots_0 = augknt([0, quantile(temp1, (1:l0-1)/l0), max(time)],k0);
    knots_q = augknt([0, quantile(temp2, (1:lq-1)/lq), max(temp2)],kq);
end

q_0 = size(knots_0, 2) - k0; 
q_q = size(knots_q, 2) - kq;

beta = est(2:p).';
theta = est(1+p:p+q_q).';
alpha = est(1+p+q_q:end).';
[~, grad_beta] = objective_func_beta(beta, x, time, delta, theta, alpha, ...
    knots_0, knots_q, k0, kq, true);
[~, grad_sieve] = objective_func_sieve([theta; alpha], x, time, delta, beta, ...
    knots_0, knots_q, k0, kq, true);
grad = [grad_beta grad_sieve];

% after training, it is possible all cumulative hazard are smaller than
% one knot, then we will get zero columns in grad.
grad(:, ~any(grad, 1)) = [];
% approximate fisher information
fish = inv((grad.' * grad));
end

