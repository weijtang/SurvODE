function main(N, seed, data_setting, knots_setting, ci)
% estimation and inference for the general linear transformation model
%   knots_setting: "K1", "K2", "K3", or "K4", one of four combinations. 
%       See mle.m for details.
%   ci: true or false - return se if true
disp(seed)


% load data
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(data_setting), '.mat');
if isfile(data_file)
    load(data_file, 'x', 'time', 'delta');
else
    generator(N, seed, data_setting);
    load(data_file, 'x', 'time', 'delta');
end
[~, p] = size(x);

tic
est = mle(x, time, delta, knots_setting);
runtime = toc;
est_r = est(1:end-1);
succ_ind = est(end);

if succ_ind
    if ci
        fish = inference(N, seed, data_setting, knots_setting, est_r.');
        se_beta = fish(1:p-1, 1:p-1);
        se_beta = sqrt(diag(se_beta));
        save(strcat('../res/survode_cgd/res_ltm_N', num2str(N), '_seed', num2str(seed), ...
            '_setting', num2str(data_setting),...
            '_knots', knots_setting, '_fish.mat'), 'se_beta')
    end
    
    res_i = [est_r; runtime].';
    save(strcat('../res/survode_cgd/res_ltm_N', num2str(N), '_seed', num2str(seed), ...
        '_setting', num2str(data_setting), '_knots', knots_setting, '.mat'), 'res_i')

end

end

