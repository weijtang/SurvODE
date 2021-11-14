function est = mle(x, time, delta, knots_setting)
%coordinate gradient descent algorithm for maximum likelihood estimation
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   knots_setting: "K1", "K2", "K3", or "K4", one of four combinations. 
%   Returns (mle, succeed indicator)

% desgin model
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

q_0 = size(knots_0, 2) - k0; 
q_q = size(knots_q, 2) - kq;

% use coxphfit as good initialization for beta
r0 = zeros(p+q_q+q_0, 1);
r0(1:p) = beta / beta(1); 

% optimize the objective function by fmincon
% history.beta = [];
history.theta = [];
history.fval = [];

options_beta = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations', 500);
options_theta = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient', true, ...
    'OutputFcn', @outfun_sieve);


% initialization
beta = r0(2:p);
theta = r0(p+1:p+q_q);
alpha = r0(p+q_q+1:end);
disp(beta)

A = [];
b = [];
lb = [];
ub = [];
nonlcon = [];
Aeq_q = zeros(1,q_q+q_0);
Aeq_q(1, q_q+1:end) = spcol(knots_0, k0, 1.5); % alpha(1.5) = 1
beq_q = log(1);

try
    for i = 1:100    
        history.theta = [];
        history.fval = [];
        % update theta and alpha
        sieve = [theta; alpha];
        fun_theta = @(r)objective_func_sieve(r, x, time, delta, beta, ...
            knots_0, knots_q, k0, kq, false);
        est_sieve = fmincon(fun_theta, sieve, A, b, Aeq_q, beq_q, lb, ub, nonlcon, ...
            options_theta);
        steptol_sieve = max(abs(est_sieve - sieve));
        theta = est_sieve(1:q_q);
        alpha = est_sieve(q_q+1:end);

        % update beta 
        fun_beta = @(r)objective_func_beta(r, x, time, delta, theta, alpha, ...
            knots_0, knots_q, k0, kq, false);
        est_beta= fminunc(fun_beta, beta, options_beta);
        steptol_beta = max(abs(est_beta - beta));
        beta =  est_beta;

        disp(beta)
        if max([steptol_beta steptol_sieve]) < 1e-3
            break
        end
    end
    succ_ind = 1;
catch
    beta = r0(2:p);
    theta = r0(p+1:p+q_q);
    alpha = r0(p+q_q+1:end);
    succ_ind = 0;
end

est = [1; beta; theta; alpha; succ_ind];



function stop=outfun_sieve(x, optimValues, ~)
stop = false;
if optimValues.iteration > 0
    steptol = max(abs(history.theta(end, :) - x.')./ (1+abs(history.theta(end, :))));
    functiontol = abs(history.fval(end) - optimValues.fval) / (1+ abs(history.fval(end)));
    constrainttol = max(optimValues.constrviolation);
    if optimValues.funccount > 100
        stop = true;
        disp('Terminate due to funccount');
    end
    if (steptol < 1e-3) && (functiontol < 1e-3) && (constrainttol < 1e-6)
        stop = true;
        disp('Terminate due to steptol and functiontol and constrainttol');
    end
end
history.fval = [history.fval; optimValues.fval];
history.theta = [history.theta; x.'];

end

end

