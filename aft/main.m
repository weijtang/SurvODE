function main(N, seed, knots_setting)
% load data
setting = 3;
if nargin < 3
    knots_setting = "quantile";
end

data_file = strcat('../data/survode/simudata_N', num2str(N), '_seed', num2str(seed), ...
    '_setting', num2str(setting), '.mat');
if isfile(data_file)
    load(data_file, 'x', 'time', 'delta');
else
    generator(N, seed, setting);
    load(data_file, 'x', 'time', 'delta');
end

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


% algorithm setting
forward = true; % if true, use forward method to calculate the gradients, otherwise use adjoint method along with parallel computing
r0 = zeros(p+q, 1);
r0(1:p) = beta; 
[est_r, runtime] = mle(x, time, delta, knots, k, r0, forward);
disp(runtime);

res_i = [est_r; runtime].';
save(strcat('../res/aft/res_aft_N', num2str(N),...
            '_seed', num2str(seed), ...
            '_knots', knots_setting, '.mat'), 'res_i')

end