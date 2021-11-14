Ns = [1 2 4 8]*1000;
m = 1000;
p = 3;
knots_setting = "quantile";
tol = 100;
px = linspace(0.01, 1.5, tol).';
resl = p+1+2*tol+2;
res_ode = zeros(length(Ns), 6);
idx = 1;
for N = Ns
    disp(N)
    l = ceil(N^(1/5));
    q = l+3;
    res = zeros(m, resl);
    se_res = zeros(m, p);
    infer = true;
    res_file = strcat('../res/cox/res_cox_summary_N', num2str(N),...
        '_knots', knots_setting,'.mat');
    if isfile(res_file)
        load(res_file, 'res', 'se_res')
    else
        for seed = 1:m
            disp(seed)
            main(N, seed, knots_setting);
            load(strcat('../res/cox/res_cox_N', num2str(N),...
                '_seed', num2str(seed), ...
                '_knots', knots_setting,'.mat'), 'res_i')
            est_r = res_i(1:p+q);
            ff = visual(seed, N, knots_setting, est_r, px, false);
            res(seed, :) = [est_r(1:p), res_i(end), ff];
            if infer
                fish = inference(N, seed, knots_setting, est_r.');
                est_se = sqrt(diag(fish(1:p, 1:p)));
                se_res(seed, :) = est_se;
            end
        end
        save(res_file, 'res', 'se_res')
    end
    est_beta = res(:,1:p);
    runtime = res(:,p+1);
    ii = p+1;
%     est_h0 = res(:,ii+1:ii+tol);
%     est_H0 = res(:,ii+tol+1:ii+2*tol);
    IMSE_h0 = res(:,ii+2*tol+1);
    IMSE_H0 = res(:,ii+2*tol+2);
    true_beta = [1, 1, 1];
    disp("beta bias:")
    disp(mean(est_beta) - true_beta);
    disp("beta SE:")
    disp(std(est_beta))
    disp('ESE:')
    disp(mean(se_res))
    up = est_beta + 1.96 * se_res;
    low = est_beta - 1.96 * se_res;
    true_para = [1 1 1];
    cp = (true_para < up) .* (true_para > low);
    disp('CP:')
    disp(mean(cp))
    if N == 1000
        runtime_base = mean(runtime);
    end
    runtime = runtime / runtime_base;
    res_ode(idx,:) = [mean(IMSE_h0), std(IMSE_h0), mean(IMSE_H0), std(IMSE_H0), mean(runtime), std(runtime)];
    idx = idx+1;
end

disp(res_ode)
res_flexsurv = csvread('../res/baseline/bdd_cox_flexsurv_all.csv');

fontsize = 12;
subplot(1,2,1);
errorbar(Ns, res_ode(:,3), res_ode(:,4),'LineWidth',4);
hold on 
errorbar(Ns, res_flexsurv(:,4), res_flexsurv(:,5), 'r', 'LineWidth',2);
hold off
xlim([900 8500])
xticks(Ns)
xticklabels({'1', '2', '4', '8'})
xlabel('Sample size (\times 10^3)')
ylabel('IMSE(\Lambda_0(\cdot))')
legend({'ODE-Cox', 'Flexsurv'}, 'Location', 'northeast')
set(gca,'FontSize',fontsize)

subplot(1,2,2);
errorbar(Ns, res_ode(:,5), res_ode(:,6),'LineWidth',2);
hold on 
errorbar(Ns, res_flexsurv(:,6), res_flexsurv(:,7), 'r', 'LineWidth',2);
hold on
x = 1000:8000;
y = x/1000;
plot(x,y,'LineStyle', '--', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth',2);
hold on
x = 1000:8000;
y = x.^2/1000^2;
plot(x,y,'LineStyle', ':', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',2);
hold off
xlim([900 8500])
xticks(Ns)
xticklabels({'1', '2', '4', '8'})
xlabel('Sample size (\times 10^3)')
ylabel('Relative computing time')
legend({'ODE-Cox', 'Flexsurv', 'Reference O(N)', 'Reference O(N^2)'}, 'Location', 'northwest')
set(gca,'FontSize',fontsize)
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gcf, 'Position', [0, 0, 600, 250]);

print(gcf,strcat('../res/fig/cox_ode_summary', '.png'),'-dpng','-r500')



