N = 8000;
tol = 100;
px = linspace(0, 2, tol).';
p = 4;
l = ceil(N^(1/5));
q = l+3;
resl = p+q+q+1+2*tol+2;
res = zeros(1000, resl);
se_res = zeros(1000, p);
for seed = 1:1000
    load(strcat('../res/cox_time_varying/res_bdd_cox_tv_N', num2str(N),...
        '_knots', num2str(l),...
        '_seed', num2str(seed), '.mat'), 'res_i');
    est = res_i(1:(end-1));
    ff = visual(seed, N, est, px, false);
    res(seed, :) = [res_i, ff];
    % inference: coverage prob
    cp_file = strcat('../res/cox_time_varying/cp_bdd_cox_tv_N', num2str(N), '.mat');
    if isfile(cp_file)
        load(cp_file,'se_res');
    else
        disp(seed)
        fish_i = inference(N, seed, est.');
        fish_beta = sqrt(diag(fish_i(1:p,1:p)));
        se_res(seed, :) = fish_beta;
    end
end
save(cp_file, 'se_res');

est_beta = res(:,1:p);
runtime = res(:,p+q+q+1);
ii = p+q+q+1;
est_h0 = res(:,ii+1:ii+tol);
est_tv = res(:,ii+tol+1:ii+2*tol);
IMSE_h0 = res(:,ii+2*tol+1);
IMSE_tv = res(:,ii+2*tol+2);


% estimation
true_beta = [1,-1,-1,1];
disp("beta bias:")
disp(mean(est_beta)-true_beta);
disp("beta se:")
disp(std(est_beta))
disp('ESE:')
disp(mean(se_res))
up = est_beta + 1.96 * se_res;
low = est_beta - 1.96 * se_res;
cp = (true_beta < up) .* (true_beta > low);
disp('CP:')
disp(mean(cp))
disp("mean IMSE h0:")
disp(mean(IMSE_h0))
disp("std IMSE h0:")
disp(std(IMSE_h0))
disp("mean IMSE tv:")
disp(mean(IMSE_tv))
disp("std IMSE tv:")
disp(std(IMSE_tv))


% compare significance between imse
est_tv_mple = csvread(strcat('../res/baseline/bdd_cox_time_varying_est_tv_N', num2str(N), '.csv'));
true_tv = sin(3*pi*px/4).';
IMSE_tv_mple = zeros(1000, 1);
for i=1:1000
    IMSE_tv_mple(i) = sum(abs(true_tv - est_tv_mple(i,:)).^2) * 2/99;
end
[h, p] = ttest(IMSE_tv, IMSE_tv_mple);
disp(p)


% visualization
fontsize = 15;
true_h0 = 0.5 * ones(size(px,1),1);
true_tv = sin(3*pi*px/4);
h0_low = quantile(est_h0, 0.025);
h0_up = quantile(est_h0, 0.975);
tv_low = quantile(est_tv, 0.025);
tv_up = quantile(est_tv, 0.975);

subplot(1,3,1);
plot(px, true_h0, 'LineWidth',4)
ylim([0 1])
hold on 
plot(px, mean(est_h0), 'r', 'LineWidth',2)
hold on 
plot(px, h0_low, 'LineStyle', '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1)
hold on 
plot(px, h0_up, 'LineStyle', '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1)
hold off
xlabel('t')
ylabel('\alpha(t)')
set(gca,'FontSize',fontsize)
legend({'True value', 'Mean estimate', '95% Percentile'})

subplot(1,3,2);
plot(px, true_tv, 'LineWidth',4)
ylim([-2 2])
hold on 
plot(px, mean(est_tv),'r', 'LineWidth',2)
hold on 
plot(px, tv_low, 'LineStyle', '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1)
hold on 
plot(px, tv_up, 'LineStyle', '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1)
hold off
xlabel('t')
ylabel('\eta(t)')
legend({'True value', 'Mean estimate', '95% Percentile'})
set(gca,'FontSize',fontsize)
set(gcf, 'Position', [0, 0, 700, 350]);
% saveas(gcf,strcat('../res/fig/cox_tv_N8000_visual', '.eps'), 'epsc')