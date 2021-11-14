function  res = visual(seed, N, est_r, px, show)
% load data from cox model 
m = csvread(strcat('../data/cox_time_varying_model/bdd_cox_N', num2str(N), '_', num2str(seed),'.csv'));
time = m(:,1);
delta = m(:,2);
% z = m(:,4);
x = m(:, [3 5 6 7]);
p = size(x,2);

% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces scales with the sample size
temp = time(delta > 0);
knots = augknt([0, quantile(temp, (1:l-1)/l), max(time)],k);
q = size(knots, 2) - k;

est_theta = est_r(p+1:p+q).';
est_alpha = est_r(p+q+1:p+q+q).';

% evaluation: visualization
B = spcol(knots, k, px);
est_h0 = exp(B * est_theta);
est_tv = B * est_alpha;
%solve ode of cum_hazard func
% y0 = 0;
% tspan = px;
% [~, cum_baseline_hazard] = ode45(@(t, y) calculate_cum_hazard(t, est_theta, ...
%     knots, k), tspan, y0);
% est_H0 = cum_baseline_hazard;
% true_H0 = 0.5 * px;
% IMSE_H0 = mean(abs(true_H0 - est_H0).^2) * 3;


% evaluation: IMSE
true_h0 = 0.5 * ones(size(px,1),1);
true_tv = sin(3*pi*px/4);
IMSE_h0 = sum(abs(true_h0 - est_h0).^2) * 2/99;
IMSE_tv = sum(abs(true_tv - est_tv).^2) * 2/99;

if show
    % visualization of baseline hazard
    plot(px, true_h0)
    hold on 
    plot(px, est_h0);
    hold off
    xlabel('t')
    ylabel('\alpha(t)')

    % visualization of time varying effects
    plot(px, true_tv)
    hold on 
    plot(px, est_tv);
    hold off
    xlabel('t')
    ylabel('\eta(t)')
    res = [];
else
    res = [est_h0.', est_tv.', IMSE_h0, IMSE_tv];
end

end

% function base_hazard = calculate_cum_hazard(t, theta, knots, k)
% B = spcol(knots, k, t);% bspines matrix
% base_hazard = exp(B * theta);
% end