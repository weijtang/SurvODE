function [l, grad] = objective_func_beta(r, x, time, delta, theta, alpha, ...
    knots_0, knots_q, k0, kq, ci)
% negative log full likelihood of the general linear transformation model
% hazard function lambda(t|x) = q(Lambda) * exp(x * beta) * alpha(t)
%   r: beta - (p-1, 1) coefficient of predictors
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   theta: (q_q, 1) current coefficient of spline bases for the function q()
%   alpha: (q_0, 1) current coefficient of spline bases for the function alpha()
%   knots_0: locations of knots for alpha()
%   knots_q: locations of knots for q()
%   k_0: splines order for alpha()
%   k_q: splines order for q()
%   ci: true or false - return gradients for computing se if true
%   Returns (loss function, gradients w.r.t. r)

[N, p] = size(x); 
beta = [1; r]; % the first element is constrained to be 1 for identifiability

%%%%%%%%  calculate loss function
multi_coef = exp(x * beta);

%solve time-transform
u = unique(time);
[~,bin] = histc(time,u);
y0 = 0;
tspan = 0;
tspan(2:length(u)+1) = u;
[~, int_alpha] = ode45(@(t, y) time_transform_func(t, alpha, ...
    knots_0, k0), tspan, y0);
int_alpha = int_alpha(2:end, :);
int_alpha = int_alpha(bin,:);
time_transform = int_alpha .* multi_coef;

%solve ode of cum_hazard function
u = unique(time_transform);
[~,bin] = histc(time_transform,u);
y0 = 0;
tspan = 0;
tspan(2:length(u)+1) = u;
[~, cum_hazard] = ode45(@(t, y) hazard_ode_func(y, theta, ...
    knots_q, kq), tspan, y0);
cum_hazard = cum_hazard(2:end, :);
cum_hazard = cum_hazard(bin,:); % sort back

% for non-decreasing of input of spcol func
u = unique(cum_hazard);
[~,bin] = histc(cum_hazard,u);

colmat = spcol(knots_q, kq, brk2knt(u, 2));
pre_Bq = colmat(1:2:end, :);
pre_dBq = colmat(2:2:end, :);
Bq = pre_Bq(bin, :); % sort back
dBq = pre_dBq(bin, :);

% for tied events and non-decreasing, unique input of spcol func
u = unique(time); 
[~, bin] = histc(time, u);
pre_B0 = spcol(knots_0, k0, u);
B0 = pre_B0(bin, :);

l1 = - (Bq * theta + x * beta + B0 * alpha).' * delta;
l2 = sum(cum_hazard);
l = (l1 + l2) / N;

%%%%%%%% calculate gradient w.r.t. beta 
x = x(:,2:p);
dd = exp(Bq * theta) .* time_transform .* x;
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

