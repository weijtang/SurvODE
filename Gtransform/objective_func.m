function [l, grad] = objective_func(r, x, time, delta, q, knots, k, rho1, r1, ci)
% negative log full likelihood under the G-transform model (rho1, r1) 
% hazard function: Lambda(t|x, beta) = G(Lambda_0(t) * exp(x * beta))
%   r: (p+q, 1) parameters (beta, theta)
%       beta: (p, 1) coefficient of predictors
%       theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   knots: locations of knots
%   k: splines order
%   rho1, r1: parameters associated with the G-transform function
%   ci: true or false - return gradients for computing se if true
%   Returns (loss function, gradients w.r.t. r)

[N, p] = size(x); 

beta = r(1:p);
theta = r(1+p:p+q);

%%%%%%%%  calculate loss function
multi_coef = exp(x * beta);

%solve ode of cum_hazard func
y0 = 0;
tspan = 1e-12;
tspan(2:N+1) = time;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, theta, ...
    knots, k), tspan, y0);
cum_baseline_hazard = cum_baseline_hazard(2:end, :);

temp = cum_baseline_hazard .* multi_coef;
[cumhaz, dcumhaz] = Gtransform(temp, rho1, r1);

% for tied events and non-decreasing, unique input of spcol func
u = unique(time); 
[~, bin] = histc(time, u);
pre_B0 = spcol(knots, k, u);
B = pre_B0(bin, :);

l1 = - (log(dcumhaz) + B * theta + x * beta).' * delta;
l2 = sum(cumhaz);
l = (l1 + l2) / N;

%%%%%%%% calculate gradient 
y0 = zeros(q, 1);
tspan = 1e-8;
tspan(2: N+1) = time;
[~ ,dtheta] = ode45(@(t, y) baseline_hazard_grad_func(t, theta, ...
    knots, k).', tspan, y0);
dtheta = dtheta(2:end, :);

if rho1 > 0
    grad_coef = - (rho1-1)./(rho1*cumhaz+1) .* delta + 1;
elseif rho1==0
    grad_coef = r1 * delta + 1;
end
grad_coef = grad_coef .* dcumhaz;

if ci
    grad = - x .* delta + x .* (grad_coef .* temp);
    grad(:, 1+p:p+q) = - B .* delta + dtheta .* (grad_coef .* multi_coef);
else
    grad = - x.' * delta + x.' * (grad_coef .* temp);
    grad(1+p:p+q) = - B.' * delta + dtheta.' * (grad_coef .* multi_coef);
    grad = grad / N;
end

end

