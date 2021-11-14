import math
import datetime
import multiprocessing as mp
import numpy as np
from lifelines import CoxPHFitter
import pandas as pd
import timeit
import scipy
from scipy.special import expit
import scipy.io

import numpy as np
from lifelines import CoxPHFitter
import pandas as pd
import timeit
import scipy
from scipy.special import expit
import scipy.io


def SPR_coef(V, delta, Z, c = 0.001, c_std = 1e-5, B = 100):
    """
    V: Duration time
    delta: Event (0: censor) 
    Z: Covariate Matrix
    c: the parameter needed for the smooth funcion: sigmoid(x/cn^(-0/5))
    The output is the point estimate of theta and the estimated standard errors of theta 
    Details could check https://academic.oup.com/biostatistics/article/8/2/197/229993 

    """

    # Fit a Cox model to set the initial values
    df = pd.DataFrame(np.column_stack((V,delta,Z)))
    df.rename(columns=dict(zip(df.columns[:2], ['Time','Event'])),inplace=True)
    cph = CoxPHFitter()
    cph.fit(df, duration_col='Time', event_col='Event')
    coef = np.array(cph.params_)
    theta_ini = (coef/coef[0])[1:] 

    n = V.shape[0]
    pairs = n*(n-1)
    k = Z.shape[1] - 1
    V_aug = V.reshape(-1,1)
    V_diff = V_aug[:, None] - V_aug[:, None].T 
    V_diff = V_diff.reshape(n,n)
    ind1 = (V_diff >= 0)*1.0
    ind2 = (V_diff <= 0)*1.0
    Z_diff1 = Z[:, :, None] - Z[:, :, None].T
    Z_diff2 = Z_diff1.T
    delta_aug1 = np.broadcast_to(delta, (n,n))
    delta_aug2 = delta_aug1.T


    def loss_theta(theta, c):
        """
        Input: initial value of theta, parameter c
        Output: (loss, gradient)
        """
        beta = np.r_[1,theta]
        beta_aug = beta.reshape(1,beta.shape[0],1)
        sigma_n = n**(-0.5)*c
        par = np.einsum('ijk,tjl-> ik', Z_diff1, beta_aug)/sigma_n
        sigmoid_mat1 = expit(par)
        sig_mat = sigmoid_mat1*(1-sigmoid_mat1)
        loss_mat = ind1 * delta_aug1 * sigmoid_mat1 
        np.fill_diagonal(loss_mat, 0)
        loss = -np.sum(loss_mat)/pairs
        grad = np.zeros(k)
        for i in range(k):
            grad_mat = delta_aug1*ind1*sig_mat*Z_diff1[:,(i+1),:]/sigma_n
            np.fill_diagonal(grad_mat, 0)
            grad[i] = -np.sum(grad_mat)/pairs
        return (loss, grad)


    def bs(theta, c, bs_seed):
        beta = np.r_[1,theta]
        beta_aug = beta.reshape(1,beta.shape[0],1)
        sigma_n = n**(-0.5)*c
        par = np.einsum('ijk,tjl-> ik', Z_diff1, beta_aug)/sigma_n
        sigmoid_mat1 = expit(par)
        sig_mat = sigmoid_mat1*(1-sigmoid_mat1)
        np.random.seed(2*bs_seed)
        W = 10*np.random.beta(0.125, 0.125, size=n)
        W_aug = W.reshape(-1,1)
        W_sum = W_aug[:, None] + W_aug[:, None].T
        W_sum = W_sum.reshape(n,n)
        loss_mat = W_sum * ind1 * delta_aug1 * sigmoid_mat1 
        np.fill_diagonal(loss_mat, 0)
        loss = -np.sum(loss_mat)/pairs
        grad = np.zeros(k)
        for i in range(k):
            grad_mat = W_sum * delta_aug1*ind1*sig_mat*Z_diff1[:,(i+1),:]/sigma_n
            np.fill_diagonal(grad_mat, 0)
            grad[i] = -np.sum(grad_mat)/pairs
        return (loss, grad) 



    time1 = timeit.default_timer()


    theta_result = scipy.optimize.minimize(loss_theta, theta_ini,  jac = True, args = (c), method = 'BFGS')
    theta_point_estimate = theta_result.x


    time2 = timeit.default_timer()

    
    bs_theta = np.zeros((B,2))
    true_vec = np.ones(2)

    for i in range(100):
        print(i)
        theta_result_bs = scipy.optimize.minimize(bs, theta_point_estimate,  jac = True, args = (c, i), method = 'BFGS')
        print(theta_result_bs.x)
        bs_theta[i,:] = theta_result_bs.x
    success_indicator = np.all(np.abs(bs_theta - true_vec) <= 3, axis = 1)
    bs_theta_success = bs_theta * success_indicator.reshape(-1,1)
    bs_theta_success = bs_theta_success[~np.all(bs_theta_success == 0, axis=1)]
    se_estimate = np.std(bs_theta_success, axis = 0)

    bs_success_rate = bs_theta_success.shape[0]/B

    time3 = timeit.default_timer()



    def tau_grad_B(theta,c):
        # Input: estimate of theta
        # Output: the n-dim vector tau, each element corresponds to each data sample
        beta = np.r_[1,theta]
        beta_aug = beta.reshape(1,beta.shape[0],1)
        sigma_n = n**(-0.5)*c
        par = np.einsum('ijk,tjl-> ik', Z_diff1, beta_aug)/sigma_n
        sigmoid_mat1 = expit(par)
        sigmoid_mat2 = sigmoid_mat1.T
        mat1 = ind1 * delta_aug1 * sigmoid_mat1 * (1 - sigmoid_mat1)
        mat2 = ind2 * delta_aug2 * sigmoid_mat2 * (1 - sigmoid_mat2)
        grad = np.zeros((n,n,k))
        for i in range(k):
            grad[:,:,i] = (mat1*Z_diff1[:,(i+1),:] + mat2*Z_diff2[:,(i+1),:])/(sigma_n)
        grad_tau = np.mean(grad, axis = 0)
        return grad_tau

    temp_B_mat1 = tau_grad_B(theta_point_estimate,c)
    temp_B_mat2 = np.einsum('ij,ik->ijk',temp_B_mat1,temp_B_mat1)
    B_mat = np.mean(temp_B_mat2, axis = 0)

    epsilon = c_std*n**(-1/4)
    u1 = np.eye(k)[:,0]
    u2 = np.eye(k)[:,1]
    grad1 = tau_grad_B(theta_point_estimate + epsilon*u1 ,c)
    grad2 = tau_grad_B(theta_point_estimate + epsilon*u2 ,c)
    A_mat = np.zeros((2, 2))
    A_mat[0,0] = np.mean((grad1[:,0] - temp_B_mat1[:,0])/epsilon)
    A_mat[1,0] = np.mean((grad1[:,1] - temp_B_mat1[:,1])/epsilon)
    A_mat[0,1] = A_mat[1,0]
    A_mat[1,1] = np.mean((grad2[:,1] - temp_B_mat1[:,1])/epsilon)
    A_mat = 0.5*A_mat
    A_mat_inv = np.linalg.inv(A_mat)
    cov_mat = A_mat_inv @ B_mat @ A_mat_inv.T
    sd_estimate = np.diag(cov_mat)**0.5
    se_estimate_sandwich = sd_estimate/(n**0.5)

    result = dict()
    result['theta'] = theta_point_estimate
    result['se'] = se_estimate
    result['se_estimate_sandwich'] = se_estimate_sandwich
    result['bs_success_rate'] = bs_theta_success.shape[0]/B
    result['theta time'] = time2 - time1
    result['se time'] = time3 - time2
    result['bs theta'] = bs_theta

    return result


# To adjust c, change the corresponding value inside this function
def calculate(seed, setting, N):
    print(seed)
    data = scipy.io.loadmat('./data/simudata_N'+str(N)+'_seed'+str(seed)+'_setting'+str(setting)+'.mat')
    temp = SPR_coef(data['time'].reshape(-1), data['delta'].reshape(-1), data['x'], c = 0.01, c_std = 1e-5, B = 100)
    return temp


if __name__ == '__main__':

    start_t = datetime.datetime.now()
    N = 1000
    setting = 1
    threshold = 1.5

    num_cores = 12
    pool = mp.Pool(num_cores)

    results = [pool.apply_async(calculate, args=(seed,setting,N)) for seed in range(1, 1001)]
    true = np.array([1,1])
    theta_estimate = np.zeros((1000, 2))
    se_estimate = np.zeros((1000, 2))
    se_estimate_sandwich = np.zeros((1000, 2))
    time = np.zeros((1000,2))
    bs_success_rate = np.zeros(1000)
    bs_theta_mat = np.zeros((1000,100,2))
    i = 0
    for p in results:
        temp = p.get()
        theta_estimate[i,:] = temp['theta']
        se_estimate[i,:] = temp['se']
        se_estimate_sandwich[i,:] = temp['se_estimate_sandwich']
        time[i,0] = temp['theta time']
        time[i,1] = temp['se time']
        bs_success_rate[i] = temp['bs_success_rate']
        bs_theta_mat[i,:,:] = temp['bs theta']
        i = i + 1


    np.save('theta_N'+str(N)+'setting'+str(setting)+'finalc001.npy',theta_estimate)
    np.save('se_N'+str(N)+'setting'+str(setting)+'finalc001.npy',se_estimate)
    np.save('se_sandwich_N'+str(N)+'setting'+str(setting)+'finalc001.npy',se_estimate_sandwich)
    np.save('time_N'+str(N)+'setting'+str(setting)+'finalc001.npy',time)
    np.save('bs_success_N'+str(N)+'setting'+str(setting)+'finalc001.npy',bs_success_rate)
    np.save('bs_theta_mat'+str(N)+'setting'+str(setting)+'finalc001.npy',bs_theta_mat)


    end_t = datetime.datetime.now()

    
    elapsed_sec = (end_t - start_t).total_seconds()
    print("total seconds: " + "{:.2f}".format(elapsed_sec) + " s")