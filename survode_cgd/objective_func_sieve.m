function [l, grad] = objective_func_sieve(r, x, time, delta, beta, ...
    knots_0, knots_q, k0, kq, ci)
% negative log full likelihood of the general linear transformation model
% hazard function lambda(t|x, beta) = q(Lambda) * exp(x * beta) * alpha(t)
%   r: (q_q+q_0, 1) parameters (theta, alpha)
%       theta: (q_q, 1) parameters in function q()
%       alpha: (q_0, 1) parameters in function alpha()
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   beta: (p, 1) current parameters for predictors
%   knots_0: locations of knots for alpha()
%   knots_q: locations of knots for q()
%   k_0: splines order for alpha()
%   k_q: splines order for q()
%   ci: true or false - return gradients for computing se if true
%   Returns (loss function, gradients w.r.t. r)

[N, ~] = size(x); 
q_0 = size(knots_0, 2) - k0;
q_q = size(knots_q, 2) - kq;

beta = [1; beta]; % the first element is constrained to be 1 for identifiability
theta = r(1:q_q);
alpha = r(q_q+1:end);

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

%solve ode of cum_hazard function and the gradients w.r.t. theta
u = unique(time_transform);
[~,bin] = histc(time_transform,u);
y0 = zeros(q_q+1,1);
tspan = 0;
tspan(2:length(u)+1) = u;
[~, res] = ode45(@(t, y) forward_odesystem_func(y, theta, ...
    knots_q, kq), tspan, y0);
res = res(2:end, :);
cum_hazard = res(bin,1); % sort back
dd_theta = res(bin,2:end);

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

%%%%%%%% calculate gradient w.r.t. alpha
y0 = zeros(q_0, 1);
tspan = 1e-8;
tspan(2: length(u)+1) = u;
[~ ,int_dalpha] = ode45(@(t, y) time_transform_grad_func(t, alpha, ...
    knots_0, k0).', tspan, y0);
int_dalpha = int_dalpha(2:end, :);
int_dalpha = int_dalpha(bin,:);
dd_alpha = exp(Bq * theta) .* multi_coef .* int_dalpha;


dd = [dd_theta dd_alpha];
grad_coef = - (dBq * theta) .* delta + 1;

if ci
    grad = - Bq .* delta;
    grad(:, q_q+1:q_q+q_0) = - B0 .* delta;
    grad = grad + dd .* grad_coef;
else
    grad = - Bq.' * delta;
    grad(q_q+1:q_q+q_0) = - B0.' * delta;
    grad = grad + dd.' * grad_coef;
    grad = grad / N;
end

end

