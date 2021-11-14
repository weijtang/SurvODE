function  res = visual(seed, N, knots_setting, est_r, px, show)
% load data
setting = 2;
if nargin < 6
    knots_setting = "quantile";
end
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
load(data_file, 'x', 'time', 'delta');
p = size(x,2);

% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces scales with the sample size
temp = time(delta > 0);
if knots_setting == "quantile"
    knots = augknt([0, quantile(temp, (1:l-1)/l), max(time)],k);
elseif knots_setting == "equal"
    knots = augknt(linspace(0,max(time),l+1), k);
end
q = size(knots, 2) - k;

est_theta = est_r(p+1:p+q).';

% evaluation: visualization
B = spcol(knots, k, px);
est_h0 = exp(B * est_theta);
% solve ode of cum_hazard func
y0 = 0;
tspan = px;
[~, cum_baseline_hazard] = ode45(@(t, y) baseline_hazard_func(t, est_theta, ...
    knots, k), tspan, y0);
est_H0 = cum_baseline_hazard;


% evaluation: IMSE
true_h0 = 2;
true_H0 = 2*px;
IMSE_h0 = sum(abs(true_h0 - est_h0).^2) * 1.5/99;
IMSE_H0 = sum(abs(true_H0 - est_H0).^2) * 1.5/99;

if show
    % visualization of baseline hazard
    plot(px, true_h0)
    hold on 
    plot(px, est_h0);
    hold off
    xlabel('t')
    ylabel('\alpha(t)')
else
    res = [est_h0.', est_H0.', IMSE_h0, IMSE_H0];
end

end
