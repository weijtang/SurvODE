function main(seed, N)
disp(seed)

% load data from cox model 
m = csvread(strcat('../data/cox_time_varying_model/bdd_cox_N', num2str(N), '_', num2str(seed),'.csv'));
time = m(:,1);
delta = m(:,2);
z = m(:,4);
x = m(:, [3 5 6 7]);
p = size(x,2);


% based on time, create knots to fit B-splines
k = 4; % order 4 cubic spline
l = ceil(size(unique(time),1)^(1/5)); % l polynomial pieces scales with the sample size
temp = time(delta > 0);
knots = augknt([0, quantile(temp, (1:l-1)/l), max(time)],k);
q = size(knots, 2) - k;


% optimize the objective function by fmincon
options = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations', 500);

fun = @(r)objective_func(r, x, z, time, delta, knots, k);
r0 = zeros(p+q+q, 1);
tic
est_r = fminunc(fun, r0, options);
runtime = toc;

res_i = [est_r; runtime].';
save(strcat('../res/cox_time_varying/res_bdd_cox_tv_N', num2str(N),...
    '_knots', num2str(l),...
    '_seed', num2str(seed), '.mat'), 'res_i')
end

