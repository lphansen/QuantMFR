---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

+++ {"id": "f5e7b1d0"}

<a href="https://colab.research.google.com/github/lphansen/RiskUncertaintyValue/blob/main/expansion_suite.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

+++

# Numerical Solution: Expansion Suite

+++ {"id": "ef624dcb"}

We previously covered the small-noise expansion method and the recursive utility framework. This section covers the usage of `Expansion Suite`, an open source Python toolbox which uses these methods to solve discrete nonlinear DSGE models.



<!-- For this notebook 

- Section 1 provides a five-minute guide to solving the DSGE model using the expansion suite, as well as how to 
    - compute shock elasticities and IRF
    - approximate and simulate variables based on model solutions. 


- Section 2 provides examples of `LinQuadVar` computations, including
    - constructions, operations, and evaluations
    - adjusting time periods
    - computing expectations -->

+++

<!-- **Table of Contents**
* [1. Five-minute guide to solving the DSGE model using the Expansion Suite](#s1)
    * [1.1 Equilibrium Condition Function](#s1-1)
    * [1.2 Model Steady State Function](#s1-2)
    * [1.3 Model Variable Dimension Array](#s1-3)
    * [1.4 Model Parameter Calibration Tuple](#s1-4)
    * [1.5 Log Consumption Growth Function](#s1-5)
    * [1.6 Solve the BY model using the Expansion Suite](#s1-6)
    * [1.7 Model Solutions](#s1-7)
    * [1.8 Approximate Variables](#s1-8)
    * [1.9 Model Simulations](#s1-9)
    * [1.10 Shock Elasticities and IRF](#s1-10)
* [2. Examples to use `LinQuadVar` in Computation](#s2)
    * [2.1 Define a `LinQuadVar` ](#s2-1)
    * [2.2 `LinQuadVar` Operations](#s2-2)
    * [2.3 `LinQuadVar` Split and Concat](#s2-3)
    * [2.4 Evaluate a `LinQuadVar`](#s2-4)
    * [2.5 Adjust Conditional Information Time](#s2-5)
    * [2.6 Compute the Expectation of time $t+1$](#s2-6) -->


```{contents}
:depth: 1
```

+++

Consider the framework previously introduced:

$$
\begin{align}
\mathbb{E}[{N_{t+1}^*Q_{t+1}^* \psi_1(X_t, J_t, X_{t+1}, J_{t+1})}|{\mathfrak{A}_t}] - \psi_2(X_t, J_t)  &= 0 \tag{1}\\
  \phi\left(X_t, X_{t+1}, J_t, {\sf q}W_{t+1}, {\sf q} \right) &= 0 
\tag{2} \end{align}
$$

The `Expansion Suite` uses the function `uncertain_expansion` to approximate a solution to the above system locally. It takes the following inputs:

```{list-table}
:header-rows: 1

* - Input
  - Type
  - Description
* - `eq`
  - *callable*
  - Function which takes state and jump variables and outputs error of equilibrium conditions, i.e. RHS of equation (1)
* - `ss`
  - *callable*
  - Function which takes model parameters and outputs steady-state values of state and jump variables
* - `var_shape`
  - *tuple of int*
  - Triple of integers $(n_J,n_X,n_W)$ where $n_J$ is the number of jump variables, $n_X$ is the number of state variables, and $n_W$ is the number of shocks.
* - `args`
  - *tuple of float or ndarray*
  - Array of preference and model parameters
* - `gc`
  - *callable*
  - Function which takes in variables and outputs log-growth of consumption
* - `args`
  - tuple of float/ndarray
  - Preference and model parameters
* - `gc`
  - *callable*
  - Function which takes in variables and outputs log-growth of consumption
* - `approach`
  - *string*
  - If equal to '1', then solve the model using approach 1 (default, see Appendix C); If equal to '2', then use approach 2 (see Appendix D). The function raises an exception if other values are provided
* - `init_util`
  - *dict, or None*
  - Dictionary for initialization of the expansion coefficients $\mu^0, \Upsilon_0^2, \Upsilon_1^2,$ and $\Upsilon_2^2$. The keys are: `mu_0`, `Upsilon_0`, `Upsilon_1`, `Upsilon_2`. If *None*, zero matrices are used
* - `iter_tol`
  - *float*
  - Iteration stops when the maximum absolute difference between the approximated variables in the current and previous iteration is below `iter_tol`
* - `max_iter`
  - *int*
  - Maximum number of iterations
```

The output is of class `ModelSolution`, which stores the coefficients for the linear-quadratic approximation for the jump and state variables; continuation values; consumption growth; and log change of measure, as well as the steady-state values of each variables. 

We will now walk through the computation using a simple example.
<br>
<br>

+++

# 1 AK Model with Recursive Utility

+++

## 1.1 Model Setup <a class="anchor" id="s1-1"></a>

Consider the following recursive utility maximization problem with AK technology:

$$
\max_{C_t, I_t} V_t = \left[ (1 - \beta) \left(C_t\right)^{1-\rho}  + \beta  \left( R_t \right)^{1-\rho} \right]^{\frac 1 {1-\rho}}\tag{3}
$$
such that

$$
C_t + I_t = \alpha K_t
$$ 

$$
\log K_{t+1} - \log K_t = \frac 1 \psi \log \left[ 1 + \psi \frac {I_t}{K_t} \right] + \beta_k Z_t - \eta_k + \sigma_k W_{t+1} 
$$ 

$$
Z_{t+1} = \exp( - \beta_z) Z_t + \sigma_z W_{t+1} 
$$ 

where $\psi$ is adjustment cost parameter which determines the efficiency of converting investment to capital, and $\eta_k$ is depreciation.
$Z_t$ governs the (scalar) conditional mean of (stochastic) technology growth while $\sigma_k$ governs the stochastic volatility of technology growth. 
$W_{t+1}$ is a two-dimensional Wiener process which provides the exogenous shocks to this model.

For the purposes of computation, we divide the first two equations by $K_t$.

+++

## 1.2 First-Order Conditions <a class="anchor" id="s1-2"></a> 

Let ${\widehat G}_t = {\widehat K}_t,$ $X_t = Z_t,$ and $D_t = \frac {I_t }{K_t}.$ We can rewrite the previous equations as: 

$$
{\widehat C}_t = \log \left( \alpha - D_t  \right) + {\widehat G}_t 
$$ 

$$
{\widehat G}_{t + 1} - {\widehat G}_t =  \frac 1 \psi \log \left[ 1 + \psi D_t \right] + \beta_k X_t - \eta_k + \sigma_k W_{t+1}
$$ 

Apply the general form of FOC and costate equations and notice in this special case $\psi^x_d = 0$, we have 

The FOC is  

$$
    \mathbb{E} \left[N^*_{t + 1}Q^*_{t + 1} \beta  \left(\frac{C_{t + 1}}{G_{t + 1}}\right)^{-\rho} \left(\frac{G_{t + 1}}{G_{t}}\right)^{1-\rho} \left(\frac{C_{t}}{G_{t}}\right)^{\rho} \frac{1}{1 + \psi D_t} \widetilde{MG}_{t+1} \right] - 1 = 0
$$ 

The law of motion of costate $\widetilde{MG}$ is 

$$
    \mathbb{E} \left[N^*_{t + 1}Q^*_{t + 1} \beta \left(\frac{C_{t + 1}}{G_{t + 1}}\right)^{-\rho} \left(\frac{G_{t + 1}}{G_{t}}\right)^{1-\rho} \left(\frac{C_{t}}{G_{t}}\right)^{\rho} \widetilde{MG}_{t+1} \right ] - \widetilde{MG}_t + \frac{C_t}{G_t} = 0
$$

Finally, we can express our system of equations as

$$
\psi_1 = 
\begin{bmatrix}
0\\ 
\beta  \left(\frac{C_{t + 1}}{G_{t + 1}}\right)^{-\rho} \left(\frac{G_{t + 1}}{G_{t}}\right)^{1-\rho} \left(\frac{C_{t}}{G_{t}}\right)^{\rho} \frac{1}{1 + \psi D_t} \widetilde{MG}_{t+1} \\
\beta \left(\frac{C_{t + 1}}{G_{t + 1}}\right)^{-\rho} \left(\frac{G_{t + 1}}{G_{t}}\right)^{1-\rho} \left(\frac{C_{t}}{G_{t}}\right)^{\rho} \widetilde{MG}_{t+1}
\end{bmatrix}
$$


$$
\psi_2 = 
\begin{bmatrix}
\frac{C_t}{G_t} + D_t - {\mathbf a}\\
1\\
\widetilde{MG}_t - \frac{C_t}{G_t}
\end{bmatrix}
$$


$$\phi=\begin{bmatrix} 
&\frac{G_{t+1}} {G_t}  -  \left( 1 + \psi D_t \right)^{\frac{1}{\psi}} + exp( \beta_k Z_t - \eta_k + \sigma_k W_{t+1} ) \\
&Z_{t+1} - \exp( - \beta_z) Z_t - \sigma_z W_{t+1} 
\end{bmatrix}$$


The model has thre
e jump variables (In the expansion suite, $J_t$ is used to represent the array of jump variables):
- $\log C_t - \log G_t$ : Log consumption over capital
- $\log D_t$ : Log investment over capital
- $\log MG_t$ : Log co-state variable

and two state variables (In the expansion suite, $X_t$ is used to represent state variables):
- $\log G_{t} - \log G_{t-1}$ : Log capital growth process
- $Z_{t}$ : Exogenous technology

+++

<br>
<br>
<br>

# 2 Inputs <a class="anchor" id="s2"></a>

## 2.1 Libraries <a class="anchor" id="s2-3"></a>

Begin by installing the following libraries and downloading `RiskUncertaintyValue`, which contains the functions required to solve the model:

```{code-cell} ipython3
import os
import sys
workdir = os.getcwd()
!git clone https://github.com/lphansen/RiskUncertaintyValue 
workdir = os.getcwd() + '/RiskUncertaintyValue'             
sys.path.insert(0, workdir+'/src')                        
import numpy as np
import autograd.numpy as anp
from scipy import optimize
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=200)
from IPython.display import display, HTML
from BY_example_sol import disp
display(HTML("<style>.container { width:97% !important; }</style>"))
import warnings
warnings.filterwarnings("ignore")

from lin_quad import LinQuadVar
from lin_quad_util import E, concat, next_period, cal_E_ww, lq_sum, N_tilde_measure, E_N_tp1
from uncertain_expansion import uncertain_expansion
np.set_printoptions(suppress=True)
```

## 2.2 Equilibrium Condition Function <a class="anchor" id="s2-2"></a>

The first input for `uncertain_expansion` requires us to define the equilibrium conditions for the model. Therefore we define a function `eq` which takes in a list of jump and state variables as inputs, and outputs the equilibrium equations in the form of equations 1 and 2. The inputs required:

```{list-table}
:header-rows: 1

* - Variable
  - Type
  - Description
  - Corresponds to
* - `Var_t`
  - array-like
  - Vector of jump and state variables in the current period
  - $(J_t, X_t)$
* - `Var_tp1`
  - array-like
  - Vector of jump and state variables in the next period
  - $(J_{t+1}, X_{t+1})$
* - `W_tp1`
  - array-like
  - Vector of shock variables in the next period
  - $(W_{t+1})$
* - `q`
  - float
  - Perturbation parameter
  - $q$
* - `mode`
  - string
  - By default, this function returns the LHS of equation (1) and (2). Set mode to 'psi1', 'psi2', or 'phi' to return the corresponding function in the equilibrium conditions
  -  
* - `args`
  - tuple of float/ndarray
  - Preference and model parameters
  -  
```

+++ {"id": "ba75562f"}


While writing out `eq`, ensure that: 

- The first components of `Var_t` and `Var_tp1` are fixed as `q_t` or `q_tp1`.
- State variables follow jump variables.
- The number of jump variables (except `q_t` and `q_tp1`) equals the number of forward-looking equilibrium conditions. 
- The number of state variables equals to the number of state evolution equations.
- All values are expressed as type `float`. Otherwise, this may cause type errors.

Note that the log change of measure is skipped in the specification of equilibrium conditions. 

Notice that the first 3 parameters and $q$ in variables should be defined and called in the given order. 

For the example model, we write:

```{code-cell} ipython3
:id: 6ae5d2d3
:tags: [hide-cell]

def eq_onecap_2d(Var_t, Var_tp1, W_tp1, q, mode, *args):
    
    # Unpack model parameters
    γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z, sigma_k, sigma_z1, sigma_z2 = args
    w1_tp1, w2_tp1 = W_tp1.ravel()

    # Unpack model variables
    q_t, cmk_t, imk_t, mk_t, Z_t, gk_t = Var_t.ravel()
    q_tp1, cmk_tp1, imk_tp1, mk_tp1, Z_tp1, gk_tp1 = Var_tp1.ravel()

    # Intermedite varibles that facilitates computation
    sdf_ex = anp.log(β) - ρ*(cmk_tp1+gk_tp1-cmk_t)
    g_dep = -α_k + beta_k*Z_t + sigma_k*w1_tp1 - 0.5 * (sigma_k**2)

    ## Forward-looking conditions
    psi1_1 = 0.
    psi1_2 = anp.exp(sdf_ex + gk_tp1) * 1/(1 + ϕ_2*anp.exp(imk_t)) * mk_tp1
    psi1_3 = anp.exp(sdf_ex + gk_tp1) * mk_tp1

    psi2_1 = -a + anp.exp(cmk_t) + anp.exp(imk_t)
    psi2_2 = 1.
    psi2_3 = mk_t - anp.exp(cmk_t)

    # State evoluion processes
    phi_1 = Z_tp1 - beta_z*Z_t - sigma_z1*w1_tp1 - sigma_z2*w2_tp1
    phi_2 = gk_tp1 - ϕ_1 * anp.log(1.+ϕ_2*anp.exp(imk_t)) - g_dep
    
    if mode == 'psi1':
        return np.array([psi1_1, psi1_2, psi1_3])
    
    return anp.array([
        psi1_1 * anp.exp(q_tp1) - psi2_1,
        psi1_2 * anp.exp(q_tp1) - psi2_2,
        psi1_3 * anp.exp(q_tp1) - psi2_3,
        phi_1, phi_2])
```

+++ {"id": "d554abc5"}

## 2.3 Steady State Function <a class="anchor" id="s2-3"></a>

This function takes model parameters as input and outputs the steady states for each variable following the variable order specified in `Var_t` and `Var_tp1`. To find the steady state, we remove the time subscripts from each variable and solve the system of equations. 

Note that the first entry of the output array will always be 0, since the steady state of `q_t` is 0.

```{code-cell} ipython3
:id: e0d871c4
:tags: [hide-cell]

def ss_onecap_2d(*args):

    γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z, sigma_k, sigma_z1, sigma_z2 = args

    def f(imk):
        g_k = ϕ_1 * np.log(1.+ ϕ_2 * np.exp(imk)) - α_k - 0.5 * (sigma_k**2)
        sdf_ex = np.log(β) - ρ*g_k
        cmk = np.log(a - np.exp(imk))
        mk = 1./(anp.exp(sdf_ex  + g_k)* 1/(1 + ϕ_2*anp.exp(imk)))
        return mk - anp.exp(cmk) - anp.exp(sdf_ex + g_k) * mk

    imk = optimize.bisect(f,-10,-2.5)
    cmk = np.log(a - np.exp(imk))
    g_k = ϕ_1 * np.log(1. + ϕ_2 * np.exp(imk)) - α_k - 0.5 * (sigma_k**2)
    sdf_ex = np.log(β) - ρ*g_k
    mk = 1./(anp.exp(sdf_ex + g_k)* 1/(1 + ϕ_2*anp.exp(imk)))

    return np.array([0., cmk, imk, mk, 0., g_k])
```

+++ {"id": "aadc21b2"}

## 2.4 Specifying Number of Variables <a class="anchor" id="s2-4"></a>

We need to specify the number of jump variables, state and shock variables as an array `(n_J, n_X, n_W)`.

```{code-cell} ipython3
:id: 6f4f50b9
:tags: [hide-cell]

var_shape = (3, 2, 2)
```

+++ {"id": "715b17d5"}

## 2.5 Model Parameters <a class="anchor" id="s2-5"></a>

Next, we need to specify all the model parameters using a Python tuple. 

- The first three parameters should be γ, β, ρ, respectively. The package is designed for recursive utlity.
- Please use 1.00001 for γ and ρ as approximation for 1.0.
- We suggest not using matrix when specifying the equilibrium conditions and parameters.

<!-- Here, we use the specification described in <insert reference here> -->

```{code-cell} ipython3
:id: 43a4ac46
:tags: [hide-cell]

σ_k = 0.477*2
sigmaz1 = 0.011*2
sigmaz2 = 0.025*2

delta = 0.01
alpha = 0.0922

a = alpha
ϕ_1 = 1/8
ϕ_2 = 8
beta_k = 0.01
beta_z1 = np.exp(-0.056)
beta_z2 = np.exp(-0.145)
beta_z3 = np.exp(-0.1)

sigma_k = σ_k*0.01
sigma_z1 = sigmaz1
sigma_z2 = sigmaz2

α_k = 0.04
β = np.exp(-delta)
γ = 12.0
ρ = 1.5 

args = (γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z1, sigma_k, sigma_z1, sigma_z2)
```

+++ {"id": "238e1c2a"}

## 2.6 Log Consumption Growth  <a class="anchor" id="s2-6"></a>

Finally, we need to specify how to compute the log-growth of consumption $\log{C_{t+1}/C_t}$ using the defined variables `Var_t` and `Var_tp1` in the equilibrium conditions.  

In this example, we simply use the decomposition:

$$
\log{C_{t+1}/C_t} = \log{C_{t+1}/K_{t+1}} + \log{K_{t+1}/K_t} - \log{C_t/K_t}
$$

- The reason to specify log consumtion growth is that we need to approximate recursive utility as described in the notes section 3.2.2 and section 3.2.3
- In habit formation models, $C_t$ is replaced with $U_t$, in which case we would specify the log growth process for $U_t$ correspondingly.

```{code-cell} ipython3
:tags: [hide-cell]

def gc_onecap_2d(Var_t, Var_tp1, W_tp1, q, *args):
    
    # Unpack model parameters
    γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z, sigma_k, sigma_z1, sigma_z2 = args

    # Unpack model variables
    q_t, cmk_t, imk_t, mk_t, Z_t, gk_t = Var_t.ravel()
    q_tp1, cmk_tp1, imk_tp1, mk_tp1, Z_tp1, gk_tp1 = Var_tp1.ravel()

    # Compute log consumption growth
    gc_tp1 = cmk_tp1 + gk_tp1 - cmk_t
    
    return gc_tp1
```

+++ {"id": "020b982b"}

## 2.7 Computing ModelSol <a class="anchor" id="s2-7"></a>

Now we are able to use steps 2.1 to 2.5 as inputs to the `uncertain_expansion` function. The solution is stored in a class named `ModelSolution`(Please refer to uncertainexpansion.ipynb for details.) under the Linear Quadratic Framework. 

See [src/uncertain_expansion.py](https://github.com/lphansen/RiskUncertaintyValue/blob/main/src/uncertain_expansion.py) for the source code of the expansion suite.

```{code-cell} ipython3
:id: 217c9448
:outputId: aab57b59-6dce-4452-dfde-a6ee003cdbd7
:tags: [hide-cell]

# Solve the One-Capital Adjustment Model
eq = eq_onecap_2d
ss = ss_onecap_2d
gc = gc_onecap_2d
var_shape = (3, 2, 2)
args = (γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z1, sigma_k, sigma_z1, sigma_z2)
approach = '1' # See Exploring Recursive Utility Appendix for difference about approach '1' and '2'
ModelSol = uncertain_expansion(eq, ss, var_shape, args, gc, approach = '1')
```

+++ {"id": "961bdc55"}

<br>
<br>
<br>

# 3 Outputs

## 3.1 List of Outputs <a class="anchor" id="s3-1"></a>

We now examine the contents of `ModelSol`, which contains the attributes listed below. Each approximation is stored in a class `LinQuadVar`, which contains the coefficients for $X_{1,t}, X_{2,t}, X_{1,t}'X_{1,t}, W_{t+1}, W_{t+1}'W_{t+1}, X_{1,t}'W_{t+1}$ and the constant.


```{list-table}
:header-rows: 1

* - Input
  - Type
  - Description
* - `JXn_t`
  - *LinQuadVar*
  - Approximation of jump and state variables at time $t$. Replace `n` with `0,1,2` to get the zeroth, first and second-order contribution. Omit `n` to get the full approximation. The variables are indexed in the order specified in Section 2. 
* - `Jn_t`
  - *LinQuadVar*
  - Same as `JXn_t` but limited to jump variables.
* - `Xn_tp1`
  - *LinQuadVar*
  - Same as `JXn_t` but limited to state variables.
* - `JXn_t_tilde`
  - *LinQuadVar*
  - Same as `JXn_t` but using distorted measure. This variation is also available for `Jn_t` and `Xn_tp1`.
* - `util_sol`
  - *dict*
  - Dictionary containing solutions of the continuation values, including $\mu^0, \Upsilon_0^2, \Upsilon_1^2,$ and $\Upsilon_2^2$ etc.
* - `vmrn_tp1`
  - *LinQuadVar*
  - Approximation of continuation values $\widehat{V}^1_{t+1}-\widehat{R}^1_t$ . Replace `n` with `0,1,2` to get the zeroth, first and second-order contribution. Omit `n` to get the full approximation. 
* - `gcn_tp1`
  - *LinQuadVar*
  - Approximation of consumption growth $\widehat{C}_{t+1}-\widehat{C}_t$ . Replace `n` with `0,1,2` to get the zeroth, first and second-order contribution. Omit `n` to get the full approximation. 
* - `ss`
  - *dict*
  - Steady states for state and jump variables
* - `log_N_tilde`
  - *LinQuadVar*
  - Approximation for the log change of measure
```

For example, we can obtain the coefficients for the first-order contribution of $\log{C_t/K_t}$ in the following way, noting that `cmk_t` was listed as the first jump variable when we specified the equilibrum conditions.

```{code-cell} ipython3
:id: 3fa57ffa
:outputId: 1cdb816b-4ed7-4c19-e009-60fdcf049004

## Log consumption over capital approximation results, shown in the LinQuadVar coefficients form
ModelSol['JX1_t'][0].coeffs
```

We can also display the full second-order approximation of $\log{C_t/K_t}$ using the `disp` function, which renders a `LinQuadVar` object in LATEX form.

```{code-cell} ipython3
## Log consumption over capital approximation results, shown in the Latex analytical form
disp(ModelSol['JX_t'][0],'\\log\\frac{C_t}{K_t}') 
```

Another example:

```{code-cell} ipython3
:id: d5b18ec2
:outputId: 1b2ebb60-d294-4281-e573-a76dcbc55b4d

## Log capital growth process second order approximation results
disp(ModelSol['X2_tp1'][0],'\\log\\frac{K_{t+1}^2}{K_t^2}') 
```

+++ {"id": "0190ee8e"}


## 3.2 Approximate Related Variables <a class="anchor" id="s3-2"></a>

We can use the method `ModelSol.approximate` to approximate variables which can be expressed as functions of the state and jump variables. For example, suppose we want to approximate log investment growth $\log{I_{t+1}/I_t}$. We can express this as:

$$
\log{I_{t+1}/I_t} = \log{I_{t+1}/K_{t+1}} + \log{K_{t+1}/K_t} - \log{I_t/K_t}
$$

We use this expression to write the function `gi_tp1_approx`, which is used as an input in `ModelSol.approximate` below.

```{code-cell} ipython3
:id: 145a47f6
:outputId: ae3e6f0d-39c8-47a6-b2b3-934eaa6cffaa
:tags: [hide-cell]

def gi_tp1_approx(Var_t, Var_tp1, W_tp1, q, *args):

    # Parameters for the model
    γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z1, sigma_k, sigma_z1, sigma_z2 = args

    # Variables: (q_t, q_tp1 is excluded when using the method in `ModelSolution`)
    cmk_t, imk_t, mkmc_t, Z_t, gk_t = Var_t.ravel()
    cmk_tp1, imk_tp1, mkmc_tp1, Z_tp1, gk_tp1 = Var_tp1.ravel()
    
    gi_tp1 = imk_tp1 + gk_tp1 -imk_t
    
    return gi_tp1

n_J, n_X, n_W = var_shape
gi_tp1_list = ModelSol.approximate(gi_tp1_approx, args)
gi_tp1 = {'gi_tp1':gi_tp1_list[0],\
        'gi0_tp1':gi_tp1_list[1],\
        'gi1_tp1':gi_tp1_list[2],\
        'gi2_tp1':gi_tp1_list[3]}
gi_tp1['gi1_tp1'].coeffs
```

+++ {"id": "bbd28ede"}

<br>

## 3.3 Simulate Variables <a class="anchor" id="s3-3"></a>
Given a series of shock processes, we can simulate the path of our state and jump variables using the `ModelSolution.simulate` method. Here, we simulate 400 periods of i.i.d standard multivariate normal shocks.

```{code-cell} ipython3
:id: 58417c59

Ws = np.random.multivariate_normal(np.zeros(n_W),np.eye(n_W),size = 400)
JX_sim = ModelSol.simulate(Ws)
```

+++ {"id": "738d5ef7"}

<br>

## 3.4 Shock Elasticities and IRF <a class="anchor" id="s3-4"></a>

`ModelSolution` contains a method `elasticities` used to compute the shock elasticities of jump and state variables. We do this by defining a `log_SDF_ex` function, which expresses the log stochastic discount factor exclusive of the change of measure ${N}^*_{t+1}$ and variable $Q_{t+1}^*$. Shock elasticities are great tools for nonlinear impulsive response analyses, and in linear models, the shock elasticities correspond to the traditional impulse response function (IRF). 

`ModelSolution` also has a method `IRF` to compute the traditional IRF of all jump and state variables. Below we calculate the IRF of all variables with respect to the first shock for 400 periods.

The figure plots elasticities and IRF of investment rate with regard to shock 1. 

```{code-cell} ipython3
:id: c2fb03dd
:tags: [hide-cell]

def log_SDF_ex(Var_t, Var_tp1, W_tp1, q, *args):

    # Log stochastic discount factor exclusive of the change of measure N and Q.

    # Parameters for the model 
    γ, β, ρ, a, ϕ_1, ϕ_2, α_k, beta_k, beta_z, sigma_k1, sigma_k2, sigma_z = args

    # Variables: (q_t, q_tp1 is excluded when using the method in `ModelSolution`) 
    cmk_t, imk_t, mk_t, Z_t, gk_t = Var_t.ravel()
    cmk_tp1, imk_tp1, mk_tp1, Z_tp1, gk_tp1 = Var_tp1.ravel()

    sdf_ex = anp.log(β) - ρ*(cmk_tp1+gk_tp1-cmk_t)
    
    return sdf_ex

expo_elas, price_elas = ModelSol.elasticities(log_SDF_ex, T=40, shock=0)
states, controls = ModelSol.IRF(T = 40, shock = 0) 
```

```{code-cell} ipython3
:tags: [hide-cell]

import matplotlib.pyplot as plt

time_periods = list(range(1, 41))  # 40 time periods starting from 1

fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # Adjust the figure size as needed

# Plotting Exposure Elasticity in the first subfigure
axs[0].plot(time_periods, expo_elas[:, 1], label='Exposure Elasticity', marker='o')
axs[0].set_xlabel('Time Period')
axs[0].set_ylabel('Exposure Elasticity')
axs[0].set_title('Exposure Elasticity Over Time')
axs[0].grid(True)

# Plotting Price Elasticity in the second subfigure
axs[1].plot(time_periods, price_elas[:, 1], label='Price Elasticity', marker='x')
axs[1].set_xlabel('Time Period')
axs[1].set_ylabel('Price Elasticity')
axs[1].set_title('Price Elasticity Over Time')
axs[1].grid(True)

# Plotting IRF in the third subfigure
axs[2].plot(time_periods, controls[:, 1], label='IRF', marker='^')
axs[2].set_xlabel('Time Period')
axs[2].set_ylabel('Impulse Response')
axs[2].set_title('IRF Over Time')
axs[2].grid(True)

plt.tight_layout()
plt.show()

```

+++ {"id": "eab36438"}

<br>
<br>

# 4 Using `LinQuadVar` in Computation <a class="anchor" id="s2"></a>

In the previous section, we saw how to use `uncertain_expansion` to approximate variables and store their coefficients as `LinQuadVar` objects. In this section, we explore how to manipulate `LinQuadVar` objects for different uses.

To aid our examples, we first extract the steady states for the state evolution processes from the previous model solution:

See [src/lin_quad.py](https://github.com/lphansen/RiskUncertaintyValue/blob/main/src/lin_quad.py) for source code of `LinQuadVar` definition.

```{code-cell} ipython3
:id: b5d7ee6b
:outputId: 882669aa-3be7-4c20-dd54-1273a7c1826d

n_J, n_X, n_W = ModelSol['var_shape']
X0_tp1 = LinQuadVar({'c':np.array([[ModelSol['ss'][1]],[ModelSol['ss'][2]]])}, shape = (2, n_X, n_W))
X0_tp1.coeffs
```

+++ {"id": "a5270666"}

## 4.1 `LinQuadVar` Operations <a class="anchor" id="s4-1"></a>
We can sum multiple LinQuads together in two different ways. Here we demonstrate this with an example by summing the zeroth, first and second order contributions of our approximation for capital growth. 

```{code-cell} ipython3
:id: 205dd148
:outputId: 75bfa88f-3b99-4535-cf0d-121591eae233

gk_tp1 = X0_tp1[0] + ModelSol['X1_tp1'][0]  + 0.5 * ModelSol['X2_tp1'][0] 
disp(gk_tp1,'\\log\\frac{K_{t+1}}{K_t}') 
```

In the next example, we sum together the contributions for both capital growth and technology:

```{code-cell} ipython3
:id: b5747b24
:outputId: b23404ea-0275-4a57-812e-a675e5f3d37d

lq_sum([X0_tp1, ModelSol['X1_tp1'], 0.5 * ModelSol['X2_tp1']]).coeffs
```

+++ {"id": "5c1952e3"}

## 4.2 `LinQuadVar` Split and Concat <a class="anchor" id="s4-2"></a>
`split` breaks multiple dimensional LinQuad into one-dimensional LinQuads, while `concat` does the inverse.

```{code-cell} ipython3
:id: d085149f
:outputId: a6bdf64b-17a6-48c5-b159-5715a34e0174

gk_tp1, Z_tp1 = ModelSol['X1_tp1'].split()
concat([gk_tp1, Z_tp1])
```

+++ {"id": "e67b04c1"}

## 4.3 Evaluate a `LinQuadVar` <a class="anchor" id="s4-3"></a>
We can evaluate a LinQuad at specific state $(X_{t},W_{t+1})$ in time. As an example, we evaluate all 5 variables under steady state with a multivariate random normal shock.

```{code-cell} ipython3
:id: 5782676a
:outputId: 0e8e4471-839b-4679-c7a2-fd8d6a8e6879

x1 = np.zeros([n_X ,1])
x2 = np.zeros([n_X ,1])
w = np.random.multivariate_normal(np.zeros(n_W),np.eye(n_W),size = 1).T
ModelSol['JX_tp1'](*(x1,x2,w))
```

+++ {"id": "5051f96c"}

## 4.4 Next period expression for `LinQuadVar` <a class="anchor" id="s4-4"></a>
`ModelSol` allows us to express a jump variable $J_t$ as a function of $t$ state and shock variables. Suppose we would like to compute its next period expression $J_{t+1}$ with shocks. The function `next_period` expresses $J_{t+1}$ in terms of $t$ state variables and $t+1$ shock variables. For example, we can express the $t+1$ expression for the first-order contribution to consumption over capital as:

```{code-cell} ipython3
:id: 7e088a47
:outputId: 116af9b7-9791-4391-b7a4-8d9ee133baff

cmk1_tp1 = next_period(ModelSol['J1_t'][0], ModelSol['X1_tp1'])
disp(cmk1_tp1, '\\log\\frac{C_{t+1}^1}{K_{t+1}^1}') 
```

```{code-cell} ipython3
:id: fa8ab664
:outputId: f4244fa2-2a8a-4a42-8a33-2ab7407c7c8d

cmk2_tp1 = next_period(ModelSol['J2_t'][0], ModelSol['X1_tp1'], ModelSol['X2_tp1'])
disp(cmk2_tp1, '\\log\\frac{C_{t+1}^2}{K_{t+1}^2}') 
```

+++ {"id": "daa0c4c1"}

## 4.6 Compute the Expectation of time $t+1$ `LinQuadVar` <a class="anchor" id="s4-6"></a>

Suppose the distribution of shocks has a constant mean and variance (not state dependent). Then, we can use the `E` function to compute the expectation of a time $t+1$ `LinQuadVar` as follows:

```{code-cell} ipython3
:id: 3140272c
:outputId: 418a75e2-ff56-4bd7-8471-67634f26654c

E_w = ModelSol['util_sol']['μ_0']
cov_w = np.eye(n_W)
E_ww = cal_E_ww(E_w, cov_w)
E_cmk2_tp1 = E(cmk2_tp1, E_w, E_ww)
disp(E_cmk2_tp1, '\mathbb{E}[\\log\\frac{C_{t+1}^2}{K_{t+1}^2}|\mathfrak{F_t}]')
```

+++ {"id": "18a1686c"}

Suppose the distribution of shock has a state-dependent mean and variance (implied by $\tilde{N}_{t+1}$ shown in the notes), we can use `E_N_tp1` and `N_tilde_measure` to compute the expectation of time $t+1$ `LinQuadVar`.

```{code-cell} ipython3
:id: a58482f7
:outputId: b8f2fc85-7c37-4d01-a5ed-c6b4d6d940a0

N_cm = N_tilde_measure(ModelSol['util_sol']['log_N_tilde'],(1,n_X,n_W))
E_cmk2_tp1_tilde = E_N_tp1(cmk2_tp1, N_cm)
disp(E_cmk2_tp1_tilde, '\mathbb{\\tilde{E}}[\\log\\frac{C_{t+1}^2}{K_{t+1}^2}|\mathfrak{F_t}]')
```
