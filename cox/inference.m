function fish = inference(N, seed, knots_setting, r)
%Construct approximate confidence interval of beta, alpha(t)
%   approximate fisher information through E(grad^2)

% load data from cox model 
setting = 1;
if nargin < 4
    knots_setting = "quantile";
end
load(strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat'), 'x', 'time', 'delta')
p = size(x,2);

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


beta = r(1:p);
theta = r(1+p:p+q);

%%%%%%%%  calculate loss function
multi_coef = exp(x * beta);

% for tied events and non-decreasing, unique input of spcol func
u = unique(time); 
[~, bin] = histc(time, u);
pre_B0 = spcol(knots, k, u);
B = pre_B0(bin, :);

%solve ode of cum_hazard func
y0 = 0;
tspan = 1e-12;
tspan(2:N+1) = time;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, theta, ...
    knots, k), tspan, y0);
cum_baseline_hazard = cum_baseline_hazard(2:end, :);

%%%%%%%% calculate gradient w.r.t. beta
grad = - x .* delta + x .* (cum_baseline_hazard .* multi_coef);
%%%%%%%% calculate gradient w.r.t. theta
dtheta1 = - B .* delta;

% solve ode of gradient w.r.t. theta
y0 = zeros(q, 1);
tspan = 1e-8;
tspan(2: N+1) = time;
[~ ,dtheta2] = ode45(@(t, y) baseline_hazard_grad_func(t, theta, ...
    knots, k).', tspan, y0);
dtheta2 = dtheta2(2:end, :);
dtheta = dtheta1 + dtheta2 .* multi_coef;

grad(:,1+p:p+q) = dtheta;

% approximate fisher information
fish = inv((grad.' * grad));
end