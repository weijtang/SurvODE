function time_transform = time_transform_func(t, alpha, knots_0, k0)
% model function alpha(t)
%   t: a scalar
%   alpha: (q_0, 1) coefficient of spline bases for the function alpha()
%   knots_0: locations of knots for alpha()
%   k_0: splines order for alpha()
B0 = spcol(knots_0, k0, t);% bspines matrix of alpha()
time_transform = exp(B0 * alpha);
end