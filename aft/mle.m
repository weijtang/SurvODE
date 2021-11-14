function [est, runtime] = mle(x, time, delta, knots, k, r0, forward)
%coordinate gradient descent algorithm for maximum likelihood estimation
%under the AFT model lambda(t|x, beta) = q(Lambda) * exp(x * beta)
%   x: (N, p) predictors
%   time: (N, 1) survival time
%   delta: (N, 1) censoring status
%   knots: locations of knots
%   k: splines order
%   r0: (p+q, 1) init parameters (beta, theta)
%       beta: (p, 1) init coefficients of predictors
%       theta: (q, 1) init coefficients of spline bases for function q
%   forward: ture or false - compute gradients via forward method if true, 
%     otherwise compute gradients via adjoint method 
%   Returns (mle, computation time)

[N, p] = size(x);

% optimize the objective function by fmincon
history.theta = [];
history.fval = [];

options_theta = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'OutputFcn', @outfun_sieve_warmup);
options_beta = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations', 500);

% initialization
beta = r0(1:p);
theta = r0(p+1:end);
disp(beta)
n_warmup = 2;

tic
for i = 1:200    
    % update theta
    if i > n_warmup
        options_theta = optimoptions('fminunc', ...
            'SpecifyObjectiveGradient', true, ...
            'OutputFcn', @outfun_sieve);
    end
    fun_theta = @(r)objective_func_sieve(r, x, time, delta, beta, ...
        knots, k, false, forward);
    est_theta = fminunc(fun_theta, theta, options_theta);
    steptol_theta = max(abs(est_theta - theta));
    theta = est_theta;
    
    % update beta 
    fun_beta = @(r)objective_func_beta(r, x, time, delta, theta, ...
        knots, k, false);
    est_beta= fminunc(fun_beta, beta, options_beta);
    steptol_beta = max(abs(est_beta - beta));
    beta = est_beta;
    
    disp(beta)
    if max([steptol_beta steptol_theta]) < 1e-3
        break
    end
end

est = [beta; theta];
runtime = toc;

function stop=outfun_sieve_warmup(x, optimValues, ~)
stop = false;
if optimValues.iteration > 0
    steptol = max(abs(history.theta(end, :) - x.')./ (1+abs(history.theta(end, :))));
    functiontol = abs(history.fval(end) - optimValues.fval) / (1+ abs(history.fval(end)));
    if optimValues.funccount > 20
        stop = true;
        disp('Terminate due to funccount');
    end
    if (steptol < 1e-2) && (functiontol < 1e-2)
        stop = true;
        disp('Terminate due to steptol and functiontol and constrainttol');
    end
end
history.fval = [history.fval; optimValues.fval];
history.theta = [history.theta; x.'];

end

function stop=outfun_sieve(x, optimValues, ~)
stop = false;
% visual(N, seed, setting, x.');

if optimValues.iteration > 0
    steptol = max(abs(history.theta(end, :) - x.')./ (1+abs(history.theta(end, :))));
    functiontol = abs(history.fval(end) - optimValues.fval) / (1+ abs(history.fval(end)));
    if optimValues.funccount > 30
        stop = true;
        disp('Terminate due to funccount');
    end
    if (steptol < max(1e-4,(100/N)^2)) && (functiontol < max(1e-4,(100/N)^2))
        stop = true;
        disp('Terminate due to steptol and functiontol and constrainttol');
    end
end
history.fval = [history.fval; optimValues.fval];
history.theta = [history.theta; x.'];

end

end

