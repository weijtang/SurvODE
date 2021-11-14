function hazard = hazard_ode_func(y, theta, knots_q, kq)
% model the ODE after the time transformation lambda(t) = q(Lambda) 
%   y: scalar
%   theta: (q_q, 1) coefficient of spline bases for the function q()
%   knots_q: locations of knots for q()
%   k_q: splines order for q()
Bq = spcol(knots_q, kq, y);% bspines matrix of q()
hazard = exp(Bq * theta);
end
