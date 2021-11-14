function  res = visual(seed, N, knots_setting, est_r, px, show)
% load data from aft model 
if nargin < 6
    knots_setting = "quantile";
end
setting = 3;
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
load(data_file, 'x', 'time', 'delta');

% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(N^(1/7));% l polynomial pieces
[beta, ~, H] = coxphfit(x, time, 'Censoring', 1-delta);
[~, bin] = histc(time, H(:,1));
bin(bin == 0) = 1;
x_coef = exp(x * beta);
temp = x_coef .* H(bin,2);
if knots_setting == "quantile"
    knots = augknt([0, quantile(temp, (1:l-1)/l), max(temp)],k);
elseif knots_setting == "equal"
    knots = augknt(linspace(0, 2 * max(temp),l+1), k);
end

p = size(x,2);
q = size(knots, 2) - k;

est_theta = est_r(p+1:p+q).';

% evaluation: visualization
B = spcol(knots, k, px);
est_q = exp(B * est_theta);


% evaluation: IMSE
true_q = 2./(1+px);
IMSE_q = sum(abs(true_q - est_q).^2) * 1.5/99;

if show
    plot(px, true_q)
    hold on 
    plot(px, est_q);
    hold off
    xlabel('t')
    ylabel('\alpha(t)')
    res = [];
else
    res = [est_q.', IMSE_q];
end

end
