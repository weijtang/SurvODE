<!-- ABOUT THE PROJECT -->

## Survival Analysis via Ordinary Differential Equation

This repo provides implementation for the proposed estimation and inference procedure in the following paper:

[Survival Analysis via Ordinary Differential Equation](https://arxiv.org/abs/2009.03449). 

[Weijing Tang](https://sites.google.com/umich.edu/weijingtang/home), [Kevin He](https://sph.umich.edu/faculty-profiles/he-zhi.html), [Gongjun Xu](https://sites.google.com/umich.edu/gongjunxu/home), and [Ji Zhu](http://dept.stat.lsa.umich.edu/~jizhu/). 

## Getting Started

### Prerequisites
The implementation of the proposed method is build on top of the following toolboxes in [Matlab](https://www.mathworks.com/products/matlab.html).  
* [Optimization Toolbox](https://www.mathworks.com/help/optim/index.html?s_tid=CRUX_lftnav) for solving constrained and unconstrained optimization problems
* [Curve Fitting Toolbox](https://www.mathworks.com/help/curvefit/index.html?s_tid=CRUX_lftnav) for constructing the sieve space using B-splines

<!-- USAGE EXAMPLES -->

## Usage

### The Cox Model

Fitting **the Cox model**: <img src="https://render.githubusercontent.com/render/math?math=\Lambda'_x(t) = \alpha(t)\exp(x^\top \beta)"> with unspecified <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)">

* The data generation: run [cox/generator.m](https://github.com/weijtang/SurvODE/blob/master/cox/generator.m) and set the setting as 1. The generated data will be saved in "data/survode/". 
* The proposed ODE-Cox approach where <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)"> is approximated by the cubic B-splines: run [cox/main.m](https://github.com/weijtang/SurvODE/blob/master/cox/main.m) to obtain the point estimate and run [cox/inference.m](https://github.com/weijtang/SurvODE/blob/master/cox/inference.m) to obtain the standard error estimate. See [cox/summary.m](https://github.com/weijtang/SurvODE/blob/master/cox/summary.m) for example.
* The comparison with the existing maximum partial likelihood estimation: the implementation [baseline/cox_mple.R](https://github.com/weijtang/SurvODE/blob/master/baseline/cox_mple.R) is based on the "coxph" function in the R package "survival".
* The comparison with the the parametric method in [Royston and Parmar (2002)](https://doi.org/10.1002/sim.1203), where the log-transformed baseline cumulative hazard is modeled as a natural cubic spline function of the log-transformed time: the implementation [baseline/cox_flexsurv.R](https://github.com/weijtang/SurvODE/blob/master/baseline/cox_flexsurv.R) is based on the "flexsurvspline" function in the R package "flexsurv".

### The Time-Varying Cox Model

Fitting **the time-varying Cox model**: <img src="https://render.githubusercontent.com/render/math?math=\Lambda'_x(t) = \alpha(t)\exp(x^\top \beta\ %2B {z}^\top \eta(t))"> with unspecified <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)"> and <img src="https://render.githubusercontent.com/render/math?math=\eta(\cdot)">

* The data generation: run [baseline/cox\_time\_varying\_generate\_data.R](https://github.com/weijtang/SurvODE/blob/master/baseline/cox_time_varying_generate_data.R) with a given sample size. The generated data will be saved in "data/cox_time_varying_model/".
* The proposed ODE-Cox-tv approach where <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)"> and <img src="https://render.githubusercontent.com/render/math?math=\eta(\cdot)"> are approximated by the cubic B-splines: see [cox\_time\_varying/main.m](https://github.com/weijtang/SurvODE/blob/master/cox_time_varying/main.m) for example.
* The comparison with the existing maximum partial likelihood estimation: the implementation [baseline/cox\_time\_varying_effects.R](https://github.com/weijtang/SurvODE/blob/master/baseline/cox_time_varying_effects.R) is based on the "coxph" function the “tt” argument set as the same cubic B-spline transformation of time in the R package "survival".

### The AFT Model

Fitting **the semi-parametric accelerated failure time (AFT) model**: <img src="https://render.githubusercontent.com/render/math?math=\Lambda'_x(t) = q(\Lambda_x(t))\exp(x^\top \beta)"> with unspecified <img src="https://render.githubusercontent.com/render/math?math=q(\cdot)">

* The data generation: run [aft/generator.m](https://github.com/weijtang/SurvODE/blob/master/aft/generator.m) and set the setting as 2. The generated data will be saved in "data/survode/". 
* The proposed ODE-AFT approach where <img src="https://render.githubusercontent.com/render/math?math=q(\cdot)"> is approximated by the cubic B-splines: see [aft/main.m](https://github.com/weijtang/SurvODE/blob/master/aft/main.m) for example. 
* I have implemented both forward and adjoint methods to compute the gradients: in [aft/main.m](https://github.com/weijtang/SurvODE/blob/master/aft/main.m), set `forward=true` for using forward method and `forward=false` for using adjoint method along with parallel computing. The default is using the forward method if the memory permits. 
* The comparison with the existing rank-based method: the implementation [baseline/aft_rank.R](https://github.com/weijtang/SurvODE/blob/master/baseline/aft_rank.R) is based on the "aftsrr" function in R package "aftgee".

### The G-Transformation Model

Fitting **the semi-parametric linear transformation model**: <img src="https://render.githubusercontent.com/render/math?math=\Lambda'_x(t) = q(\Lambda_x(t))\exp(x^\top \beta)\alpha(t)"> with specified <img src="https://render.githubusercontent.com/render/math?math=q(\cdot | \rho, r)">  and unspecified <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)">. In particular, the solution of the initial value problem <img src="https://render.githubusercontent.com/render/math?math=\tilde{\Lambda}'_x(t) = q(\tilde{\Lambda}_x(t) \big|\rho, r)"> with <img src="https://render.githubusercontent.com/render/math?math=\tilde{\Lambda}(0) = 0"> corresponds to the G-transformation function in [Zeng and Lin(2007)](https://doi.org/10.1111/j.1369-7412.2007.00606.x), which covers the class of Box-Cox transformations and the class of logarithmic transformations.  

* The data generation: run [Gtransform/generator.m](https://github.com/weijtang/SurvODE/blob/master/Gtransform/generator.m) and set the setting as 3. The generated data will be saved in "data/survode/". 
* The proposed ODE-LT approach where <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)"> is approximated by the cubic B-splines: see [Gtransform/main.m](https://github.com/weijtang/SurvODE/blob/master/Gtransform/main.m) for example.
* The comparison with the existing nonparametric maximum likelihood estimation in [Zeng and Lin(2007)](https://doi.org/10.1111/j.1369-7412.2007.00606.x): the implementation [baseline/NPMLE_Trasnform/main.m](https://github.com/weijtang/SurvODE/blob/master/baseline/NPMLE_Transform/main.m) is modified based on the original implementation using the EM algorithm in [Zeng and Lin(2007)](https://doi.org/10.1111/j.1369-7412.2007.00606.x).

### The Linear Transformation Model

Fitting **the general (or nonparametric) linear transformation model**: <img src="https://render.githubusercontent.com/render/math?math=\Lambda'_x(t) = q(\Lambda_x(t))\exp(x^\top \beta)\alpha(t)"> with both <img src="https://render.githubusercontent.com/render/math?math=\alpha(\cdot)"> and <img src="https://render.githubusercontent.com/render/math?math=q(\cdot)"> unspecified

* The data generation: run [survode_cgd/generator.m](https://github.com/weijtang/SurvODE/blob/master/survode_cgd/generator.m) and set the setting between 1 and 4. In particular, the Cox, the AFT, and the logarithmic transformation model are mis-specified under setting 4. The generated data will be saved in "data/survode/". 
* The proposed ODE-Flex approach: see [survode_cgd/main.m](https://github.com/weijtang/SurvODE/blob/master/survode_cgd/main.m) and [survode_cgd/mle.m](https://github.com/weijtang/SurvODE/blob/master/survode_cgd/mle.m) for example.
* The comparison with the existing smoothed partial rank method in [Song et al. (2006)](https://doi.org/10.1093/biostatistics/kxl001): the implementation [baseline/SPR_LT.py](https://github.com/weijtang/SurvODE/blob/master/baseline/SPR_LT.py) is able to reproduce the simulation results in Song et al. (2006). Note that SPR introduces an additional parameter c in the objective function to improve the estimation accuracy. We evaluated SPR with various values of the parameter c and the sample size N under our data settings 1)-4) in the manuscript. 

### Knots Placements

I have implemented two placements of knots for cubic splines. The argument `knots_setting` in "main.m" can be set as:

* `knots_setting = "quantile"`: the interior knots arelocated at the quantiles of the distinct observation time points.
* `knots_setting = "equal"`: the interior knots equally separate the time interval from 0 to the maximum of observed times.

The default placement is `knots_setting = "quantile"`.

<!-- CONTACT -->

## Contact

Weijing Tang (University of Michigan) - weijtang@umich.edu

