function [l, grad] = objective_func(r, x, z, time, delta, knots,k)
% negative log full likelihood under the time-varying cox model
% hazard function lambda(t|x, beta) = lambda_0(t) * exp(x * beta + z * eta(t))
%   r: (p+q+q, 1) parameters (beta, theta, alpha)
%       beta: (p, 1) coefficients of predictors
%       theta: (q, 1) coefficients of spline bases for the baseline hazard function
%       alpha: (q, 1) coefficients of spline bases for the time varying coefficient
%   x: (N, p) predictors
%   z: (N, 1) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   knots: locations of knots
%   k: splines order
%   Returns (loss function, gradients w.r.t. r)

[N, p] = size(x); 
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

l1 = - (B * theta + (B * alpha).* z + x * beta).' * delta;


%solve ode of cum_hazard function
y0 = zeros(N,1);
tspan = 1e-12;
tspan(2) = 1;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, theta, ...
    alpha, z, time, knots, k), tspan, y0);
cum_baseline_hazard = cum_baseline_hazard(end, :);
l2 = cum_baseline_hazard * multi_coef;
l = (l1 + l2) / N;

%%%%%%%% calculate gradient w.r.t. beta
grad = - x.' * delta + x.' * (cum_baseline_hazard.' .* multi_coef);

%%%%%%%% calculate gradient w.r.t. theta and alpha
d1 = - B.' * delta;
d1(q+1:q+q) = - B.' * (delta.*z);

% solve ode of gradient w.r.t. theta and alpha
y0 = zeros(q+q, 1);
tspan = 1e-8;
tspan(2) = 1;
[~ ,d2] = ode45(@(t, y) baseline_hazard_grad_func(t, theta, ...
    alpha, z, multi_coef, time, knots, k), tspan, y0);
d2 = d2(end,:);
dd = d1 + d2.';

grad(1+p:p+q+q) = dd;

grad = grad / N;
% disp("done")

end