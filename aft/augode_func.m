function res = augode_func(y, theta, x_coef, knots, k)
% return the augmented ODE for computing gradients
%   y: [cumhaz adjoint d_cumhaz_theta]
%   theta: (q, 1) coefficients of spline bases for function q
%   x_coef: exp(x * beta)
%   knots: locations of knots
%   k: splines order
colmat = spcol(knots, k, brk2knt(y(1), 2));
B = colmat(1, :);
dB = colmat(2, :);
res_1 = exp(B * theta) * x_coef;
res_2 = - y(2) * res_1 * (dB * theta);
res_4 = - y(2) * res_1 * B;
res = [res_1, res_2, res_4].';
end

