function res = visual(N, seed, setting, knots_setting, est, pxt, pxq, show)
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
load(data_file, 'x', 'time', 'delta');

[N, p] = size(x);
k0 = 4; % order 4 cubic spline
kq = 4; 
l0 = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces for alpha()
lq = ceil(N^(1/7)); % l polynomial pieces for q()

temp1 = time(delta > 0);
[beta, ~, H] = coxphfit(x, time, 'Censoring', 1-delta);
[~, bin] = histc(time, H(:,1));
bin(bin == 0) = 1;
x_coef = exp(x * beta);
temp2 = x_coef .* H(bin,2);

if knots_setting == "K1"
    knots_0 = augknt(linspace(0,max(time),l0+1), k0);
    knots_q = augknt(linspace(0, 2 * max(temp2),lq+1), kq);
elseif knots_setting == "K2"
    knots_0 = augknt(linspace(0,max(time),l0+1), k0);
    knots_q = augknt([0, quantile(temp2, (1:lq-1)/lq), max(temp2)],kq);
elseif knots_setting == "K3"
    knots_0 = augknt([0, quantile(temp1, (1:l0-1)/l0), max(time)],k0);
    knots_q = augknt(linspace(0, 2 * max(temp2),lq+1), kq);
elseif knots_setting == "K4"
    knots_0 = augknt([0, quantile(temp1, (1:l0-1)/l0), max(time)],k0);
    knots_q = augknt([0, quantile(temp2, (1:lq-1)/lq), max(temp2)],kq);
end

% visualization bound 
up = pxt(end);
qup = pxq(end);

q_0 = size(knots_0, 2) - k0;
q_q = size(knots_q, 2) - kq;
est_theta = est(p+1:p+q_q).';
est_alpha = est(p+q_q+1:p+q_q+q_0).';


% estimator of q() and alpha()
Bq = spcol(knots_q, kq, pxq);
est_q = exp(Bq * est_theta);
B0 = spcol(knots_0, k0, pxt);
est_a = exp(B0 * est_alpha);


if setting == 1 % Cox model correctly specified
    true_a = pxt.^3;
    true_q = ones(size(pxq, 1), 1);
    est_a = est_a * 1.5^3;
    est_q = est_q / 1.5^3;
    yup = up.^3+0.5;
elseif setting == 2 % G-transform model & AFT model
    true_a = ones(size(pxt, 1), 1) * 2;
    true_q = exp(-pxq);
    est_a = est_a * 2;
    est_q = est_q / 2;
    yup = 3;
elseif setting == 3 % AFT model
    true_a = ones(size(pxt, 1), 1);
    true_q = 2./(1+pxq);
    yup = 2;
elseif setting == 4 % General linear transformation model
    true_a = log(1+pxt);
    true_q = log(1+pxq) +2;
    est_a = est_a * log(2.5);
    est_q = est_q / log(2.5);
    yup = log(1+up)+0.5;
end




if show
    % visualize alpha
    subplot(1,2,1);
    plot(pxt, true_a, 'LineWidth',2)
    xlim([0 up])
    ylim([0 yup])
    hold on 
    plot(pxt, est_a, 'r', 'LineWidth',2)
    hold off
    xlabel('t')
    ylabel('\alpha(t)')
    % visualize q
    subplot(1,2,2);
    plot(pxq, true_q, 'LineWidth',2)
    xlim([0 qup])
    hold on 
    plot(pxq, est_q,'r', 'LineWidth',2)
    hold off
    xlabel('u')
    ylabel('q(u)')
    res = [];
else
    res = [est_a.', est_q.'];
end
end