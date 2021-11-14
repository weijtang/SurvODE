function [l, grad] = objective_func_sieve(r, x, time, delta, beta, ...
    knots, k, ci, forward)
% negative log full likelihood of the aft model
% hazard function lambda(t|x, beta) = q(Lambda) * exp(x * beta)
%   r: (q, 1) init coefficients of spline bases for function q
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   beta: (p, 1) current coefficients of predictors
%   knots: locations of knots
%   k: splines order
%   ci: true or false - return gradients for computing se if true
%   forward: ture or false - compute gradients via forward method if true, 
%     otherwise compute gradients via adjoint method 
%   Returns (loss function, gradients w.r.t. theta)
[N, ~] = size(x); 
q = size(knots, 2) - k;

theta = r;

%%%%%%%%  calculate loss function
x_coef = exp(x * beta);
temp = time .* x_coef;

u = unique(temp);
[~,bin] = histc(temp,u);

if forward
    %solve augmented ode of cum_hazard and gradient w.r.t. theta
    y0 = zeros(q+1,1);
    tspan = 0;
    tspan(2:length(u)+1) = u;
    [~, res] = ode45(@(t, y) forward_odesystem_func(y, theta, ...
        knots, k), tspan, y0);
    res = res(2:end, :);
    cum_hazard = res(bin,1); % sort back
    dd = res(bin,2:end);
else %adjoint
    %solve ode of cum_hazard function
    y0 = 0;
    tspan = 0;
    tspan(2:length(u)+1) = u;
    [~, cum_hazard] = ode45(@(t, y) hazard_ode_func(y, theta, ...
        knots, k), tspan, y0);
    cum_hazard = cum_hazard(2:end, :);
    cum_hazard = cum_hazard(bin,:); % sort back
    %solve augmented ode of gradient w.r.t. theta
    dd = zeros(N, q);
    parfor (i = 1:N, 200)
        y0 = [cum_hazard(i), 1, zeros(1, q)].';
        tspan = time(i);
        tspan(2) = 0;
        [~, d] = ode45(@(t, y) augode_func(y, theta, ...
            x_coef(i), knots, k), tspan, y0);
        dd(i,:) = d(end,3:end);
    end
%     y0 = [cum_hazard, ones(N,1), zeros(N,q)]';
%     tspan = [1, 0];
%     [~, d] = ode45(@(t, y) adjoint_odesystem_func(y, theta, x_coef, time, knots, k), tspan, y0);
%     dd = reshape(d(end,:), [], N).';
%     dd = dd(:, 3:end);
end

u = unique(cum_hazard);
[~,bin] = histc(cum_hazard,u);
colmat = spcol(knots, k, brk2knt(u, 2));
pre_Bq = colmat(1:2:end, :);
pre_dBq = colmat(2:2:end, :);
Bq = pre_Bq(bin, :); % sort back
dBq = pre_dBq(bin, :);
grad_coef = - (dBq * theta) .* delta + 1;

if ci
    grad = - Bq .* delta;
    grad = grad + dd .* grad_coef;
else
    grad = - Bq.' * delta;
    grad = grad + dd.' * grad_coef;
    grad = grad / N;
end

l1 = - (Bq * theta + x * beta).' * delta;
l2 = sum(cum_hazard);
l = (l1 + l2) / N;

end

