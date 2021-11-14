function generator(N, seed, setting, rho1, r1)
%generate data from linear transformation model
%   specify beta = [1, 1, 1]
%   specify q(u)
%   specify alpha(t)

%generate covariates x: one binary and two continuous
%generate uniform u

if nargin > 3
    rho1 = rho1; r1=r1;
else
    rho1=0; r1=1;
end

p = 3;
beta = [1; 1; 1];
rho = 0.5^2;
mu = 0;
rng(seed);
pd = makedist('Normal','mu',mu,'sigma',rho);
tnorm = truncate(pd,-2,2);
x = random(tnorm, N, p);
u = rand(N, 1);

if setting == 1 % Cox model correctly specified
    true_hazard_ode = @(y, t, x, beta) 1 * exp(x*beta)*t^3;
    up = 5; % censoring rate 25%
elseif setting == 2 % G-transform model & AFT model
    if rho1 > 0
        true_hazard_ode = @(y, t, x, beta) power(rho1*y+1,1-1/rho1) * exp(x*beta) * t^2 * 2;
    elseif rho1==0
        true_hazard_ode = @(y, t, x, beta) exp(- r1*y + x*beta)  * 2;
    end
    up = 4; % censoring rate 27%
elseif setting == 3 % AFT model
    true_hazard_ode = @(y, t, x, beta) 2/(1+y) * exp(x*beta);
    up = 3.5; % censoring rate 25%
elseif setting == 4 % General linear transformation model
    true_hazard_ode = @(y, t, x, beta) (log(1+y) + 2) * exp(x*beta) * log(1+t);
    up=4; % censoring rate 23%
end

pre_time = zeros(N, 1);
for i = 1:N
    y0 = 1e-8;
    tspan = [0 up+1];
    sol = ode45(@(t, y) true_hazard_ode(y, t, x(i,:), beta), tspan, y0);
    cum_hazard = @(t) deval(sol, t);
    f = @(t) cum_hazard(t) + log(1-u(i,:));
    if f(up+1) < 0
        pre_time(i,:) = up+1;
    else
        pre_time(i,:) = fzero(f, [0 up+1]);
    end
end     

%generate censor time
pre_censoring = up * rand(N, 1);
time = min(pre_time, pre_censoring);
delta = double(pre_censoring >= pre_time);

% hist(pre_time)
disp(quantile(pre_time(), [0.33 0.66]))
disp(1-mean(delta))

[time, idx] = sort(time);
delta = delta(idx);
x = x(idx, :);

save(strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat'), 'x', 'time', 'delta')


end
