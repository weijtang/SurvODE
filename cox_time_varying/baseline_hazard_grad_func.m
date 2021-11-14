function res = baseline_hazard_grad_func(t, theta, alpha, z, multi_coef, time, knots, k)
% return the gradient of lambda_0(t) * exp(z * eta(t)) w.r.t. (theta, alpha)
%   t: a number between 0 and 1
%   theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   alpha: (q, 1) coefficients of spline bases for the time varying coefficient
%   z: (N, 1) predictors
%   multi_coef: exp(x * beta)
%   time: (N, 1) survival time
%   knots: locations of knots
%   k: splines order
q = size(theta, 1);
B = spcol(knots, k, t * time);
baseline_hazard = exp(B * theta + (B * alpha) .* z);
temp = baseline_hazard .* multi_coef .* time;
baseline_hazard_grad = temp.' * B;
baseline_hazard_grad(:, q+1:q+q) = (temp .* z).' * B;
res = baseline_hazard_grad.';

end