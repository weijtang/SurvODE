function [l, grad] = objective_func(r, x, time, delta, p, q, knots,k)
% negative log full likelihood under the cox model
% hazard function lambda(t|x, beta) = lambda_0(t) * exp(x * beta)
%   r: (p+q, 1) parameters (beta, theta)
%       beta: (p, 1) coefficients of predictors
%       theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   knots: locations of knots
%   k: splines order
%   Returns (loss function, gradients w.r.t. r)

[N, ~] = size(x); 

beta = r(1:p);
theta = r(1+p:p+q);

%%%%%%%%  calculate loss function
multi_coef = exp(x * beta);

% for tied events and non-decreasing, unique input of spcol func
u = unique(time); 
[~, bin] = histc(time, u);
pre_B0 = spcol(knots, k, u);
B = pre_B0(bin, :);

l1 = - (B * theta + x * beta).' * delta;

%solve ode of cum_hazard func
y0 = 0;
tspan = 1e-12;
tspan(2:N+1) = time;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, theta, ...
    knots, k), tspan, y0);
cum_baseline_hazard = cum_baseline_hazard(2:end, :);

l2 = cum_baseline_hazard.' * multi_coef;
l = (l1 + l2) / N;

%%%%%%%% calculate gradient w.r.t. beta
grad = - x.' * delta + x.' * (cum_baseline_hazard .* multi_coef);

%%%%%%%% calculate gradient w.r.t. theta
dtheta1 = - B.' * delta;

% solve ode of gradient w.r.t. theta
y0 = zeros(q, 1);
tspan = 1e-8;
tspan(2: N+1) = time;
[~ ,dtheta2] = ode45(@(t, y) baseline_hazard_grad_func(t, theta, ...
    knots, k).', tspan, y0);
dtheta2 = dtheta2(2:end, :);
dtheta = dtheta1 + dtheta2.' * multi_coef;

grad(1+p:p+q) = dtheta;
grad = grad / N;
end

