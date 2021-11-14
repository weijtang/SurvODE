function time_transform_grad = time_transform_grad_func(t, alpha, knots_0, k0)
% return the gradient of function alpha(t) w.r.t. alpha
%   t: a scalar
%   alpha: (q_0, 1) coefficient of spline bases for the function alpha()
%   knots_0: locations of knots for alpha()
%   k_0: splines order for alpha()
B0 = spcol(knots_0, k0, brk2knt(t, 1));
time_transform_grad = exp(B0 * alpha) .* B0;
end