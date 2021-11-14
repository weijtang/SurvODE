settings = [1 2 3 4];
Ns = [1 2 4 8]*1000;
knots_settings = ["K4"];% "K3" "K2" "K1"];
for knots_setting = knots_settings
    for setting = settings
        for N = Ns
            disp(setting)
            disp(N)
            m = 1000;
            p = 3;
            resl = p+1;
            res = zeros(m, resl);
            succ_seed = false(m,1);
            se_res = zeros(m, p-1);
            res_file = strcat('../res/survode_cgd/res_survode_summary_N', num2str(N), ...
                '_setting', num2str(setting),'_knots', knots_setting,'.mat');
            if isfile(res_file)
                load(res_file, 'res', 'se_res', 'succ_seed')
            else
                parfor seed = 1:m
                    try
                        disp(seed)
                        main(N, seed, setting, knots_setting, false);
                        res_i = load(strcat('../res/survode_cgd/res_ltm_N', num2str(N),...
                                '_seed', num2str(seed),'_setting', num2str(setting),...
                                '_knots', knots_setting, '.mat'));
                        res(seed, :) = res_i.res_i([1:p end]);
                        fish = inference(N, seed, setting, knots_setting, res_i.res_i(1:(end-1)));
                        se_beta = fish(1:p-1, 1:p-1);
                        se_res(seed, :) = sqrt(diag(se_beta));
                    catch 
                        continue; 
                    end
                    succ_seed(seed) = true;
                end
                save(res_file, 'res', 'se_res', 'succ_seed')
            end
            disp(strcat('success rate:', num2str(sum(succ_seed)),'/', num2str(m)))
            res = res(succ_seed,:);
            se_res = se_res(succ_seed,:);

            est_beta = res(:,2:p);
            runtime = res(:,p+1);
            disp('runtime:')
            disp(mean(runtime))
            ii = p+1;
            IMSE_q = res(:,end);
            true_beta = [1, 1];
            disp("beta bias:")
            disp(mean(est_beta) - true_beta);
            disp("beta SE:")
            disp(std(est_beta))
            disp('ESE:')
            disp(mean(se_res))
            up = est_beta + 1.96 * se_res;
            low = est_beta - 1.96 * se_res;
            true_para = [1 1];
            cp = (true_para < up) .* (true_para > low);
            disp('CP:')
            disp(mean(cp))
        end
    end
end

settings = [1 2 3 4];
Ns = [1 2 4 8]*1000;

