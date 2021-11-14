function baseline_hazard = baseline_hazard_func(time, theta, knots, k)
% model baseline hazard function lambda_0(t)
%   time: (N, 1) survival time
%   theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   knots: locations of knots
%   k: splines order
B = spcol(knots, k, brk2knt(time, 1));
baseline_hazard = exp(B * theta);
end

