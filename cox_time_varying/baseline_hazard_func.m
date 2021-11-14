function baseline_hazard = baseline_hazard_func(t, theta, alpha, z, time, knots, k)
% model the time-varying part of the hazard function lambda_0(t) * exp(z * eta(t))
%   t: a number between 0 and 1
%   theta: (q, 1) coefficients of spline bases for the baseline hazard function
%   alpha: (q, 1) coefficients of spline bases for the time varying coefficient
%   z: (N, 1) predictors
%   time: (N, 1) survival time
%   knots: locations of knots
%   k: splines order
B = spcol(knots, k, t * time);% bspines matrix
baseline_hazard = exp(B * theta + (B * alpha).* z) .* time;
end