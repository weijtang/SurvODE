function res = adjoint_odesystem_func(y, theta, x_coef, time, knots, k)
% model the augmented ODE system for computing gradients simultaneously for
% n subjects
%   y: [cumhaz adjoint d_cumhaz_theta]
%   theta: (q, 1) coefficients of spline bases for function q
%   x_coef: exp(x * beta)
%   time: (N, 1) observed time points
%   knots: locations of knots
%   k: splines order

% change the size of y to be: number of equations by number of initial (2+q, N)
[N, ~] = size(x_coef);
y = reshape(y, [], N);

% compute the spline and its gradient
u = unique(y(1,:).'); 
[~, bin] = histc(y(1,:).', u);
pre_B = spcol(knots, k, brk2knt(u, 2));
B = pre_B(bin*2-1, :);
dB = pre_B(bin*2, :);

res_1 = exp(B * theta) .* x_coef .* time;
res_2 = - y(2,:).' .* res_1 .* (dB * theta);
res_4 = - y(2,:).' .* res_1 .* B;
res = [res_1, res_2, res_4].';
res = res(:);
end

