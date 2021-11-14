function [l, grad] = objective_func_beta(r, x, time, delta, theta, ...
    knots, k, ci)
% negative log full likelihood of aft model
% hazard function lambda(t|x, beta) = q(Lambda) * exp(x * beta)
%   r: (p, 1) init coefficients of predictors
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   theta: (q, 1) current coefficients of spline bases for function q
%   knots: locations of knots
%   k: splines order
%   ci: true or false - return gradients for computing se if true
%   Returns (loss function, gradients w.r.t. r)
[N, ~] = size(x); 

beta = r;
%%%%%%%%  calculate loss function
x_coef = exp(x * beta);
temp = time .* x_coef;

u = unique(temp);
[~,bin] = histc(temp,u);

%solve ode of cum_hazard function
y0 = 0;
tspan = 0;
tspan(2:N+1) = u;
[~, cum_hazard] = ode45(@(t, y) hazard_ode_func(y, theta, ...
    knots, k), tspan, y0);
cum_hazard = cum_hazard(2:end, :);
cum_hazard = cum_hazard(bin,:); % sort back

u = unique(cum_hazard);
[~,bin] = histc(cum_hazard,u);

colmat = spcol(knots, k, brk2knt(u, 2));
pre_Bq = colmat(1:2:end, :);
pre_dBq = colmat(2:2:end, :);
Bq = pre_Bq(bin, :); % sort back
dBq = pre_dBq(bin, :);

l1 = - (Bq * theta + x * beta).' * delta;
l2 = sum(cum_hazard);
l = (l1 + l2) / N;

%%%%%%%% calculate gradient w.r.t. beta
dd = exp(Bq * theta) .* temp .* x;
grad_coef = - (dBq * theta) .* delta + 1;

if ci
    grad = - x .* delta;
    grad = grad + dd .* grad_coef;
else
    grad = - x.' * delta;
    grad = grad + dd.' * grad_coef;
    grad = grad / N;
end


end

