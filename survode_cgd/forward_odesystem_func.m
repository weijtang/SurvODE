function res = forward_odesystem_func(y, theta, knots, k)
% model the augmented ODE for computing gradients via forward method and cumulative hazard
%   y: [cumhaz cumhaz_theta]
%   theta: (q, 1) coefficient of spline bases for the function q()
%   knots: locations of knots for q()
%   k: splines order for q()
colmat = spcol(knots, k, brk2knt(y(1), 2));
B = colmat(1, :);
dB = colmat(2, :);
d_cumhaz = exp(B * theta);
d_cumhaz_theta = d_cumhaz * (B.' + dB * theta * y(2:end));
res = [d_cumhaz; d_cumhaz_theta];
end

