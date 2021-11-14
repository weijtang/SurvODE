function main(N, seed, knots_setting)

% load data
setting = 2;
if nargin < 3
    knots_setting = "quantile";
end
rho1 = 0; r1 = 1;
data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');


if isfile(data_file)
    load(data_file, 'x', 'time', 'delta');
else
    generator(N, seed, setting, rho1, r1);
    load(data_file, 'x', 'time', 'delta');
end
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


% optimize the objective function by fmincon
options = optimoptions('fminunc', ...
    'Algorithm','trust-region', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations', 500);

fun = @(r)objective_func(r, x, time, delta, q, knots, k, rho1, r1, false);
r0 = zeros(p+q,1);
tic
est_r = fminunc(fun, r0, options);
runtime = toc;
% disp(est_r);

res_i = [est_r; runtime].';
save(strcat('../res/Gtransform/res_Gtransform_N', num2str(N),...
    '_seed', num2str(seed), ...
    '_knots', knots_setting,'.mat'), 'res_i')

end