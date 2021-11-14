Ns = [1 2 4 8]*1000;
m = 1000;
p = 3;
knots_setting = "quantile";
tol = 100;
px = linspace(0.01, 1.5, tol).';
resl = p+1+tol+1;
res_ode = zeros(length(Ns), 4);
idx = 1;
for N = Ns
    disp(N)
    l = ceil( N^(1/7));
    q = l+3;
    res = zeros(m, resl);
    succ_seed = false(m,1);
    se_res = zeros(m, p);
    infer = true;
    res_file = strcat('../res/aft/res_aft_summary_N', num2str(N),...
        '_knots', knots_setting,'.mat');
    if isfile(res_file)
        load(res_file, 'res', 'se_res', 'succ_seed')
    else
        parfor seed = 1:m
            try
                disp(seed)
                res_seed = strcat('../res/aft/res_aft_N', num2str(N),...
                    '_seed', num2str(seed), ...
                    '_knots', knots_setting, '.mat');
%                 main(N, seed, knots_setting);
%                 load(res_seed, 'res_i')
                if isfile(res_seed)
                    res_i = load(res_seed, 'res_i');
                else
                    main(N, seed, knots_setting);
                    res_i = load(res_seed, 'res_i');
                end
                est_r = res_i.res_i(1:p+q);
                ff = visual(seed, N, knots_setting, est_r, px, false);
                res(seed, :) = [est_r(1:p), res_i.res_i(end), ff];
                if infer
                    fish = inference(N, seed, knots_setting, est_r.');
                    est_se = sqrt(diag(fish(1:p, 1:p)));
                    se_res(seed, :) = est_se;
                end
            catch 
                continue;  % Jump to next iteration of: for i
            end
            
            succ_seed(seed) = true;
        end
        save(res_file, 'res', 'se_res', 'succ_seed')
    end
    disp(strcat('fail: ', num2str(m-sum(succ_seed)),'/', num2str(m)))
    res = res(succ_seed,:);
    se_res = se_res(succ_seed,:);
    
    est_beta = res(:,1:p);
    runtime = res(:,p+1);
    ii = p+1;
    IMSE_q = res(:,end);
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
    temp= sort(runtime);
    runtime = temp(51:(sum(succ_seed)-50));
    res_ode(idx,:) = [mean(IMSE_q), std(IMSE_q), mean(runtime), std(runtime)];
    idx = idx+1;
end

disp(res_ode)
res_aftgee = csvread('../res/baseline/bdd_aftgee_runtime.csv');

fontsize = 12;
subplot(1,2,1);
errorbar(Ns, res_ode(:,1), res_ode(:,2),'LineWidth',2);
hold off
xlim([900 8500])
xticks(Ns)
xticklabels({'1', '2', '4', '8'})
xlabel('Sample size (\times 10^3)')
ylabel('IMSE(q(\cdot))')
legend({'ODE-AFT'}, 'Location', 'northeast')
set(gca,'FontSize',fontsize)

subplot(1,2,2);
errorbar(Ns, res_ode(:,3), res_ode(:,4),'LineWidth',2);
hold on 
errorbar(Ns, res_aftgee(:,2), res_aftgee(:,3), 'r', 'LineWidth',2);
hold on
x = 1000:8000;
y = x/1000;
plot(x,y,'LineStyle', '--', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth',2);
hold on
x = 1000:8000;
y = x.^2/1000^2;
plot(x,y,'LineStyle', '--', 'Color', [0.9290 0.6940 0.1250], 'LineWidth',2);
hold off
xlim([900 8500])
ylim([0.2 150])
xticks(Ns)
xticklabels({'1', '2', '4', '8'})
xlabel('Sample size (\times 10^3)')
ylabel('Relative computing time')
legend({'ODE-AFT', 'Rank-based', 'Reference O(N)', 'Reference O(N^2)'}, 'Location', 'northwest')
set(gca,'FontSize',fontsize)
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gcf, 'Position', [0, 0, 600, 250]);

print(gcf,strcat('../res/fig/aft_ode_summary', '.png'),'-dpng','-r500')






