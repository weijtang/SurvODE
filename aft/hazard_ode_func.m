function hazard = hazard_ode_func(y, theta, knots, k)
% model ODE function q(u)
%   y: cumulative hazard
%   theta: (q, 1) coefficients of spline bases for function q
%   knots: locations of knots
%   k: splines order
B = spcol(knots, k, y);% bspines matrix
hazard = exp(B * theta);
end

