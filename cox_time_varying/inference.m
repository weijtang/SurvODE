function fish = inference(N, seed, r)
%Construct approximate confidence interval of r
%   approximate fisher information through E(grad^2)

% load data from cox model 
m = csvread(strcat('../data/cox_time_varying_model/bdd_cox_N', num2str(N), '_', num2str(seed),'.csv'));
time = m(:,1);
delta = m(:,2);
z = m(:,4);
x = m(:, [3 5 6 7]);
p = size(x,2);

% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces scales with the sample size
temp = time(delta > 0);
knots = augknt([0, quantile(temp, (1:l-1)/l), max(time)],k);
q = size(knots, 2) - k;

beta = r(1:p);
theta = r(1+p:p+q);
alpha = r(p+q+1:p+q+q);

%%%%%%%%  calculate loss function
multi_coef = exp(x * beta);

% for tied events and non-decreasing, unique input of spcol func
u = unique(time); 
[~, bin] = histc(time, u);
pre_B0 = spcol(knots, k, u);
B = pre_B0(bin, :);

%solve ode of cum_hazard function
y0 = zeros(N,1);
tspan = 1e-12;
tspan(2) = 1;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, theta, ...
    alpha, z, time, knots, k), tspan, y0);
cum_baseline_hazard = cum_baseline_hazard(end, :);

%%%%%%%% calculate gradient w.r.t. beta
grad = - x .* delta + x .* (cum_baseline_hazard.' .* multi_coef);

%%%%%%%% calculate gradient w.r.t. theta and alpha
z = [ones(N, 1), z];
theta = reshape(r(1+p:end),[q,size(z,2)]);
% solve ode of gradient w.r.t. theta and alpha
dd = zeros(N, q*size(z,2));
parfor (ii = 1:N, 10)
    if mod(ii, 100)==0
        disp(ii)
    end
    y0 = zeros(q*size(z,2), 1);
    tspan = 1e-8;
    tspan(2) = 1;
    [~ ,d2] = ode45(@(t, y) single_grad_func(t, theta, ...
        z(ii,:), multi_coef(ii), time(ii), knots, k), tspan, y0);
    dd(ii,:) = d2(end,:) - reshape(B(ii,:).' * z(ii,:) .* delta(ii), [1, q*size(z,2)]);
end

grad = [grad, dd];

% approximate fisher information
fish = inv((grad.' * grad));
end

function res = single_grad_func(t, theta, z, multi_coef, time, knots, k)
q = size(theta, 1);
B = spcol(knots, k, t*time);
baseline_hazard = exp(dot(B * theta, z, 2));
temp = baseline_hazard .* multi_coef .* time;
res = reshape(B.' * (temp.*z), [q*size(z,2), 1]);
end

