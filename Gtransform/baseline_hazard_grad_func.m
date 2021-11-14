function baseline_hazard_grad = baseline_hazard_grad_func(time, theta, knots, k)
% return the gradient of baseline hazard function lambda_0(t) w.r.t. theta
%   time: (N, 1) survival time
%   theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   knots: locations of knots
%   k: splines order
B = spcol(knots, k, brk2knt(time, 1));
baseline_hazard_grad = exp(B * theta) .* B;
end

