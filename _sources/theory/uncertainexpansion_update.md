---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3.9.12 ('base')
  language: python
  name: python3
---

<a href="https://colab.research.google.com/github/lphansen/RiskUncertaintyValue/blob/main/uncertainexpansion.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

+++ {"id": "4bc82505"}

# Uncertain Expansion Theory

A nonlinear DSGE model typically cannot be solved quasi-analytically.  We consider  a small-noise expansion method to find an approximate solution modified to incorporate recursive utility formulation motivated by robustness or risk concerns. 

Please refer to [*Exploring Recursive Utility*](https://larspeterhansen.org/wp-content/uploads/2023/01/2023-01-19-recursive_expansion_LPH.pdf)  note for detailed derivation and explanation. 


<br>

+++ {"id": "2420fc4c"}

<!-- # 1 Framework for Solving DSGE
A DSGE model's equilibrium conditions can be described by a system of possibly nonlinear, expectational difference equations:

$$
\begin{align}
\mathbb{E}[{N_{t+1}^*Q_{t+1}^* \psi_1(X_t, J_t, X_{t+1}, J_{t+1})}|{\mathfrak{A}_t}] - \psi_2(X_t, J_t)  &= 0 \tag{1}\\
  \phi\left(X_t, X_{t+1}, J_t, {\sf q}W_{t+1}, {\sf q} \right) &= 0 
\tag{2} \end{align}
$$
where $X_t$ is an $n$-dimensional date $t$ state vector, $J_t$ is an $m$-dimensional date $t$ jump vector, and $W_{t+1}$ is distributed as a $k$-dimensional standard normal vector. 

You can think of ${\sf q}$ as a parameter which scales the shock $W_{t+1}$. The stochastic processes $X$ and $J$ depend implicitly on ${\sf q}$. We solve for equilibrium by imposing stochastic stability, i.e. ${\sf q} = 0$ .

The first set of $m$ equations is forward looking and hence involves expectations along with contributions $N_{t+1}^*$ and $Q_{t+1}^*$ from recursive utility.  For instance, the first-order necessary condition for optimal investment might contribute a forward-looking equation with the investment to capital ratio as a jump variable.

The second set of $n$ equations gives the dynamic evolution of the state vector. For instance, the law of motion for capital contributes one such equation. Part of our desired solution is to express the jump vector as a function of the state vector.  -->


+++

# 1 Small Noise Expansion

Consider the following class of stochastics processes indexed by a scalar perturbation parameter $q$: 

$$
X_{t + 1}(q) = \psi[X_t(q), qW_{t + 1}, q] 
$$ 

We parameterize this familty so that $q = 1$ gives the model of interest. 

Denote the zero-order expansion $q = 0$ limit as: 

$$
X_{t+1}^0 = \psi (X_t^0, 0, 0)
$$

and assume there exists a second-order expansion of $X_t$ around $q = 0$: 

$$
X_t \approx X_t^0 + q X_t^1 + \frac{q^2}{2}X_t^2
$$

Processes $X_t^j,j=0,1,2$ have a recursive structure: the stochastic process $X_t^0$ can be computed first, then the process $X_t^1$ next (it depends on $X_t^0$), and then the process $X_t^2$ (it depends on both $X_t^0$ and $X_t^1$). 

$$
X_{t+1}^1 =  \begin{bmatrix} \psi_{x'}^1 \cr \psi_{x'}^2 \cr \vdots \cr \psi_{x'} ^n \end{bmatrix} X_t^1 + \begin{bmatrix} \psi_{w'}^1 \cr \psi^2_{w'}
\cr \vdots \cr \psi^n_{w'} \end{bmatrix} W_{t+1} + \begin{bmatrix} \psi^1_{\sf q}  \cr \psi^2_{\sf q}  \cr \vdots \cr \psi^n_{\sf q}  \end{bmatrix}
$$

$$
\begin{align} 
X_{t+1}^2 =  \hspace{.2cm} & \psi_{x'} X_t^2 + \begin{bmatrix} X_t^{1'} \psi_{xx'}^1 X_t^1 \cr X_t^{1'} \psi_{xx'}^2  X_t^1\cr \vdots \cr X_t^{1'} \psi_{xx'}^n X_t^1 \end{bmatrix} +  2 \begin{bmatrix} X_t^{1'} \psi_{xw'}^1 W_{t+1}  \cr X_t^{1'} \psi_{xw'}^2  W_{t+1} \cr \vdots \cr X_t^{1'} \psi_{xw'}^n W_{t+1}  \end{bmatrix} + \begin{bmatrix} {W_{t+1}}' \psi_{ww'}^1 W_{t+1}  \cr {W_{t+1}}' \psi_{ww'}^2  W_{t+1} \cr \vdots \cr {W_{t+1}}' \psi_{ww'}^n W_{t+1}  \end{bmatrix} \\ 
& +  2 \begin{bmatrix} \psi^1_{{\sf q}  x'}X_{t}^1 \cr  \psi^2_{{\sf q}  x'}X_{t}^1 \cr \vdots \cr \psi^n_{{\sf q}  x'} X_t^1 \end{bmatrix}
+ 2 \begin{bmatrix} \psi_{{\sf q}  w' }^1W_{t+1} \cr \psi_{{\sf q}  w' }^2 W_{t+1} \cr \vdots \cr \psi^n_{{\sf q}  w'} W_{t+1} \end{bmatrix}
+ \begin{bmatrix} \psi^1_{{\sf q}  {\sf q} } \cr \psi^2_{{\sf q} {\sf q} } \cr \vdots \cr \psi^n_{{\sf q} {\sf q} } \end{bmatrix}
\end{align}
$$

Let $C$ denote consumption and ${\widehat C}$ the logarithm of consumption.    Suppose that the logarithm of consumption evolves as:

$$
\begin{equation*}
{\widehat C}_{t+1}  - {\widehat C}_t = \kappa (X_t, {\sf q} W_{t+1}, {\sf q}).  
\end{equation*}
$$ 

Approximate this process by: 

$$ 
\begin{equation}
{\widehat C}_{t+1} - {\widehat C}_t \approx 
{\widehat C}_{t+1}^0 - {\widehat C}_t^0 + {\sf{q}} \left({\widehat C}_{t+1}^1 - {\widehat C}_t^1 \right)
+ \frac{{\sf q}^2} 2 \left({\widehat C}_{t+1}^2 - {\widehat C}_t^2\right) 
\end{equation}
$$

where 

$$
\begin{align*}
{\widehat C}_{t+1}^0  - {\widehat C}_t^0  = \hspace{.2cm}  & \kappa (X_t^0, 0, 0) := \eta_0^c \\ 
{\widehat C}_{t+1}^1  - {\widehat C}_t^1  = \hspace{.2cm} & \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1} + \kappa_q \\ 
{\widehat C}_{t+1}^2  - {\widehat C}_t^2  =   \hspace{.2cm}  & \kappa_{x'} X_t^2 + {X_{t}^1}' \kappa_{x,x'} X_t^1 + 2
{X_{t}^1}' \kappa_{xw'} W_{t+1} + {W_{t+1}}'\kappa_{ww'}W_{t+1} \\ 
&+ 2 \kappa_{q,x'} X_t^1 +  2 \kappa_{qw'}W_{t+1}  + \kappa_{qq} .
\end{align*}
$$

+++ {"id": "aafa5967"}

# 2 Approximating ${\widehat V}_t$ 

Consider the **recursive utility** preference specification as follows: 

$$
V_t = \left[ (1 - \beta) \left(C_t\right)^{1-\rho}  + \beta  \left( R_t \right)^{1-\rho} \right]^{\frac 1 {1-\rho}}
$$

where 

$$ 
R_t = \left( {\mathbb E} \left[ \left( V_{t+1} \right)^{1-\gamma} \mid {\mathfrak A}_t \right] \right)^{\frac 1 {1-\gamma}} .
$$ 

Note that $0<\beta<1$ is a subjective discount factor, $\rho$ is the inverse elasticity of intertemporal substitution and $\gamma$ describes risk aversion. Since continuation values are determined only up to an increasing transformation, we work with logs of the above equations:

$$
{\widehat  V}_t = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t] + \beta \exp \left[(1-\rho) {\widehat  R}_t \right] \right]
$$

where

$$
{\widehat  R}_t = {\frac 1 {1 - \gamma}} \log {\mathbb E} \left(  \exp \left[ (1 - \gamma) {\widehat V}_{t+1} \right] \mid {\mathfrak A}_t \right).
$$

While $\gamma$ is usually interpreted as a measure of risk aversion, we now turn to an alternative interpretation which allows us to incorporate ambiguity aversion into our models. (Please refer to Section 3.1 in [*Exploring Recursive Utility*](https://larspeterhansen.org/wp-content/uploads/2023/01/2023-01-19-recursive_expansion_LPH.pdf)  note)


Construct 

$$
\begin{align}
N_{t+1}^* & \stackrel{\text { def }}{=} \exp \left[(1 - \gamma) \left({\widehat V}_{t+1} - {\widehat R}_t\right) \right]  \\
Q_{t+1}^* & \stackrel{\text { def }}{=} \exp \left[(1 - \rho) \left({\widehat V}_{t+1} - {\widehat R}_t\right) \right]. 
\end{align}
$$

## 2.1 Zeroth-order contribution

We assume that the first difference in the logarithm of consumption is stationary. For models with production, this representation will become part of the derived expansion. For order zero, we presume a constant growth rate:

$$
{\widehat C}_{t+1}^0 - {\widehat  C}_t^0  = \eta_c
$$

The order-zero approximation, ${\widehat  V}_t^0 - {\widehat  C}_t^0 =\eta_{v-c}$, is determined by 

$$
\exp\left[(1-\rho) {\left( \eta_{v - c} \right) }\right] = {\frac {1 - \beta} { 1 - \beta \exp \left[ (1 - \rho) \eta_c \right]}} .
$$





## 2.2 First-order contribution

The first-order approximation of ${\widehat V}_t$ net of ${\widehat C}_t^1$ is given by 

$$
{\widehat  V}_t^1 - {\widehat  C}_t^1 = \lambda \left( {\widehat  R}_t^1 - {\widehat  C}_t^1 \right)
$$

where

$$
{\widehat  R}_t^1 - {\widehat  C}_t^1 =  \left({\frac 1 {1 - \gamma_o}} \right) \log {\mathbb  E}  \left( \exp \left[ (1 - \gamma_o) \left[ \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1}^1\right) + 
\left({\widehat C}_{t+1}^1 - {\widehat C}_{t}^1
\right)   \right] \right] \mid {\mathfrak A}_t \right)\
$$

and $0 < \lambda < 1$  satisfies 

$$
\begin{align*}
\lambda  & = \beta \exp \left[(1-\rho) \eta_c \right]
\end{align*}
$$

Combining these equations gives the following equation to be solved:

$$
\begin{align*} 
{\widehat  V}_t^1 - {\widehat  C}_t^1
= & \lambda \left({\frac 1 {1 - \gamma_o}} \right) \log {\mathbb  E}  \left( \exp \left[ (1 - \gamma_o) \left[ \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1}^1\right) + 
\left({\widehat C}_{t+1}^1 - {\widehat C}_{t}^1
\right)   \right] \right] \mid {\mathfrak A}_t \right)
 \end{align*} 
$$ 

We guess and verify that 

$$
\hat V_t^1 - \hat C _t ^1 = v_1'X_t^1 + v_0
$$

where 

$$
v_1 = \lambda(I - \lambda\psi_{x'})^{-1} \kappa_{x'}
$$

$$
v_0 = \frac{\lambda}{(1 - \lambda)}(v_1' \psi_q + \kappa_q) + \frac{\lambda(1 - \gamma_0)}{2(1 - \lambda)}|v_1'\psi_{w'} + \kappa_{w'}|^2
$$

We have 

$$
{\widehat V}_{t+1}^1 - {\widehat C}_t^1  = {\upsilon_1}' X_{t+1}^1 + \upsilon_0 + \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1} .
$$

and 

$$
{\widehat V}_{t+1}^1- {\widehat R}_t^1 = \left( {\upsilon_1}'\psi{_w'} + \kappa_{w'} \right) W_{t+1} - \left( \frac { 1 - \gamma_o } 2 \right)
\vert  {\upsilon_1}'\psi{_w'} + \kappa_{w'} \vert^2 
$$

## 2.3 Second-order contribution

The second-order approximation is: 

$$
{\widehat  V}_t^2 - {\widehat  C}_t^2 = \lambda \left( {\widehat  R}_t^2 - {\widehat  C}_t^2 \right) +  (1- \rho) (1-\lambda) \lambda  \left( {\widehat  R}_t^1 - {\widehat  C}_t^1\right)^2
$$

where

$$
{\widehat  R}_t^2  =  E \left( N_{t+1}^0  {\widehat V}^2_{t+1} \mid {\mathfrak A}_t \right),
$$

and $N_{t+1}^0$ is given by

$$
N_{t+1}^0 = \exp \left[\left(1-\gamma_o\right)\left(\widehat{V}_{t+1}^1-\widehat{R}_t^1\right)\right]
$$

Combining these formulas, gives the following equation to be solved:

$$
\begin{equation}
\begin{split}
{\widehat  V}_t^2 - {\widehat  C}_t^2 = & \lambda 
{\mathbb E}  \left[ N_{t+1}^0 \left({\widehat V}^2_{t+1} - {\widehat C}_{t+1}^2 \right) \mid {\mathfrak A}_t \right]  + \lambda  {\mathbb E}  \left[ N_{t+1}^0 \left( {\widehat C}^2_{t+1} - {\widehat C}_{t}^2 \right) \mid {\mathfrak A}_t \right] \\ 
& + \frac{(1- \rho) (1-\lambda)}{\lambda } \left[ \left( {\widehat  V}_t^1 - {\widehat  C}_t^1\right)^2 \right]  
\end{split}
\end{equation}
$$

The formula for $N_{t+1}^0$ also happens to be the zeroth-order approximation to $N_{t+1}^*$ given by definition.  Since $N_{t+1}^0$ is positive and has conditional mean equal to one, it induces a change in probability distribution whereby $W_{t+1}$ is distributed as a multivariate normal with mean $\mu^o$ and covariance matrix $I$.  The algorithm uses this change of probability distribution in the computations.  The first-order contribution to $N_{t+1}^*$ is  

$$
N_{t+1}^1=\left(\frac{1-\gamma_o}{2}\right) N_{t+1}^0\left(\widehat{V}_{t+1}^2-\widehat{R}_t^2\right)
$$



The algorithm also uses the following approximation formulas for $Q_t^*$. 

$$
\begin{aligned}
Q_{t+1}^0 &= \exp \left[(\rho-1)\left(\widehat{V}_{t+1}^0-\widehat{R}_t^0\right)\right]=1 \\
Q_{t+1}^1 &= \left.\frac{d}{d {\sf q}} \exp \left[(\rho-1)\left(\widehat{V}_{t+1}-\widehat{R}_t\right)\right]\right|_{\mathbf{q}=0}=(\rho-1)\left(\widehat{V}_{t+1}^1-\widehat{R}_t^1\right) \\
Q_{t+1}^2 &= \left.\frac{d^2}{d \mathsf{q}^2} \exp \left[(\rho-1)\left(\widehat{V}_{t+1}-\widehat{R}_t\right)\right]\right|_{\mathbf{q}=0}=(\rho-1)^2\left(\widehat{V}_{t+1}^1-\widehat{R}_t^1\right)^2+(\rho-1)\left(\widehat{V}_{t+1}^2-\widehat{R}_t^2\right)
\end{aligned}
$$

The code will produce the coefficient inputs into the following formulas:

$$
\widehat{V}_{t+1}^1-\widehat{R}_t^1=\frac{1}{1-\gamma_o}\left[\mu^0 \cdot\left(W_{t+1}-\mu^0\right)+\frac{1}{2} \mu^0 \cdot \mu^0\right]
$$

and 

$$
\begin{aligned}
\widehat{V}_{t+1}^2-\widehat{R}_t^2= & \frac{1}{2}\left(W_{t+1}-\mu^0\right)^{\prime} \Upsilon_2^2\left(W_{t+1}-\mu^0\right)-\frac{1}{2} \operatorname{tr}\left(\Upsilon_2^2\right)+\left(W_{t+1}-\mu^0\right)^{\prime}\left(\Upsilon_1^2 X_t^1+\Upsilon_0^2\right)
\end{aligned}
$$

<br><br>

+++

# 3 Solving a Planner's Problem with Recursive Utility 

Write a triangular system with stochastic growth as:

$$
\begin{align*}
X_{t+1}\left( \mathsf{q}\right) & =\psi^x \left[D_t\left( \mathsf{q}\right),  X_{t}\left( \mathsf{q}\right), {\sf q} W_{t+1}\right] 
 \cr
\log G_{t+1} \left( \mathsf{q}\right) - \log G_t \left( \mathsf{q}\right) & = \hat \psi^g\left[ D_t\left( \mathsf{q}\right), X_{t}\left( \mathsf{q}\right),  {\sf q} W_{t+1} \right] ,
\end{align*}
$$

where $D_t$ is a date $t$ decision vector for the planner.  In addition, we impose 

$$
{\widehat C}_t \left( \mathsf{q}\right)  = \hat \kappa \left[D_t \left( \mathsf{q}\right),  X_{t}\left( \mathsf{q}\right) \right] + {\widehat G}_t \left( \mathsf{q}\right) .
$$ 

In what follows for notational convenience, we will leave the ${\sf q}$ implicit.  
We use homogeneity to rewrite the utility recursion: 

$$
\begin{equation*}
{\widehat V}_t - {\widehat G}_t= {\frac 1 {1-\rho}} \log \left[ (1 - \beta) \exp\left[(1-\rho)( {\widehat C}_t - {\widehat G}_t )\right]    + \beta  \exp\left[(1-\rho)  ({\widehat R}_t -  {\widehat G}_t)\right]\right]  
\end{equation*}
$$ 

where 

$$
\begin{equation*}
{\widehat R}_t - {\widehat G}_t = \frac 1 {1-\gamma} \log \left( {\mathbb E} \left[ \exp\left( (1 - \gamma) \left[\left( {\widehat V}_{t+1}  - {\widehat G}_{t+1} \right) +  \left({\widehat G}_{t+1} - {\widehat G}_t \right)\right]  \right)  \mid {\mathfrak A}_t \right] \right) .
\end{equation*} 
$$ 

First-order conditions for $D$:

$$
\begin{align*}
& (1-\beta) \exp\left[(1- \rho ) \left({\widehat C}_t - {\widehat V}_t \right)\right] {\hat \kappa}_d( D_t, X_t) 
\cr &+ 
\beta \exp\left[(1 - \rho) \left({\widehat R}_t -  {\widehat V}_t  \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1-\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right) \right]   \left[\psi_{d'}^x ( D_t, X_t, W_{t+1} )\right]' MX_{t+1} \mid {\mathfrak A}_t \right) \cr
&+  \beta \exp\left[(1 - \rho ) \left({\widehat R}_t - {\widehat V}_t \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1-\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right) \right]   {\hat \psi}_{d'}^g ( D_t, X_t, W_{t+1} )  \mid {\mathfrak A}_t \right)  \cr
& = 0.
\end{align*}
$$ 

where $MX_{t+1}$ is the co-state, or the partial derivation of next periods value function evaluated at the next period state vector.

The corresponding co-state equation is

$$
\begin{align*}
MX_t = & (1 - \beta)  \exp\left[(1 - \rho ) \left({\widehat C}_t - {\widehat V}_t \right)\right] {\hat \kappa}_{x}(D_t,X_t)   \cr
& + \beta \exp\left[(1 - \rho ) \left({\widehat R}_t - {\widehat V}_t \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1-\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right) \right]   \left[\psi_{x'}^x ( D_t, X_t, W_{t+1} )\right]' MX_{t+1} \mid {\mathfrak A}_t \right) \cr
& +  \beta \exp\left[(1 - \rho ) \left({\widehat R}_t - {\widehat V}_t \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1-\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right) \right]   {\hat \psi}_{x}^g ( D_t, X_t, W_{t+1} ) \mid {\mathfrak A}_t \right).  
\end{align*} 
$$ 

Similarly,

$$ 
\begin{align*}
M{\widehat G}_t = & (1 - \beta)  \exp\left[(1 - \rho ) \left({\widehat C}_t - {\widehat V}_t \right)\right]    \cr
& + \beta \exp\left[(1 - \rho ) \left({\widehat R}_t - {\widehat V}_t \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1-\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right) \right] M{\widehat G}_{t+1} \mid {\mathfrak A}_t \right).  
\end{align*} 
$$ 

Notice that $M{\widehat G}_t = 1$ solves this equation.  


The solution is to find $MX_t$ and $D_t$ as time-invariant functions of $X_t,$ making the dynamic system behave in a stochastically stable manner.  

+++


To solve the system, we first transform the costate equations. 
Define the adjusted costates as: 

$$
\begin{align*}
 \widetilde{MX}_t &= \frac {MX_t  \exp\left[(1 - \rho ){\widehat V}_t + \rho \widehat{C}_t - \widehat{G_t}  \right]}  { (1 - \beta)} \\
\widetilde{MG}_t &=  \frac {M {\widehat G}_t\exp\left[(1 - \rho ){\widehat V}_t + \rho \widehat{C}_t - {\widehat G}_t \right]}  { (1 - \beta)}. 
\end{align*} 
$$ 

Then we have the modified co-state equation:

$$ 
\begin{align*}
& \widetilde{MX}_t =  \kappa_{x}(D_t,X_t)    \cr
& +
\beta{\mathbb E} 
\left(
\exp\left[(\rho-\gamma) \left({\widehat {V}}_{t+1}  - {\widehat {R}}_t \right)  -\rho (\widehat{C}_{t+1}-\widehat{C}_t) \right]  \left(\frac {G_{t+1}}{G_t} \right)\psi_{x'}^x ( D_t, X_t, W_{t+1} )' \widetilde{MX}_{t+1} \mid {\mathfrak A}_t \right) \cr
& + \beta{\mathbb E} 
\left(
\exp\left[(\rho-\gamma) \left({\widehat {V}}_{t+1}  - {\widehat {R}}_t \right)  -\rho (\widehat{C}_{t+1}-\widehat{C}_t) \right]  \psi_{x}^g ( D_t, X_t, W_{t+1} )'  \widetilde{MG}_{t+1}\mid {\mathfrak A}_t \right).  
\end{align*} 
$$ 

According to the definition of $\widetilde{MG}_t$, we have

$$
\begin{align*}
\widetilde{MG}_t = \kappa(D_t, X_t) + {\mathbb E} \left(
\exp\left[(\rho-\gamma) \left({\widehat {V}}_{t+1}  - {\widehat {R}}_t \right)  -\rho \left({\widehat {C}}_{t+1} - {\widehat {C}}_t\right)  \right] \left(\frac { G_{t+1}}{ G_t}\right)  \widetilde{MG}_{t+1}\mid {\mathfrak A}_t \right)
\end{align*}
$$ 

Write the first-order conditions as: 

$$ 
\begin{align*}
& 0 =   \kappa_d( D_t, X_t) 
\cr &+ 
\beta 
{\mathbb E} 
\left(
\exp\left[(\rho -\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right)  - \rho \left(\widehat{C}_{t+1}-\widehat{C}_t\right) \right]   \left(\frac {G_{t+1}}{G_t}\right) \psi_{d'}^x ( D_t, X_t, W_{t+1} )' {\widetilde {MX}} _{t+1} \mid {\mathfrak A}_t \right) \cr
&+  \beta 
{\mathbb E} 
\left(
\exp\left[(\rho -\gamma) \left({\widehat V}_{t+1}  - {\widehat R}_t \right)  - \rho (\widehat{C}_{t+1}-\widehat{C}_t) \right]   \psi_{d'}^g ( D_t, X_t, W_{t+1} ) {\widetilde {MG}}_{t+1}  \mid {\mathfrak A}_t \right) .
\end{align*}
$$ 

Write:

$$
\frac {S_{t+1}}{S_t} = N_{t+1}^* Q_{t+1} \beta \exp \left[ - \rho  \left( {\widehat C}_{t+1} - {\widehat C}_t \right) \right]
$$

where

$$
\begin{align*} 
N_{t+1}^* & = \exp\left[  (1-\gamma_o) \left({\widetilde V}_{t+1} - {\widetilde  R}_t\right) \right] \cr
Q_{t+1}^* & = \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right]
\end{align*}
$$ 

are terms that are contributed by recursive utility.  

Note that the FOC and costate equations can be rewritten as: 

$$
\begin{align*}
& \widetilde{MX}_t =  \kappa_{x}(D_t,X_t)    \cr
& +
\beta{\mathbb E} 
\left(
N^*_{t+1}Q^*_{t + 1}  \left(\frac {G_{t+1}}{G_t} \right)\psi_{x'}^x ( D_t, X_t, W_{t+1} )' \widetilde{MX}_{t+1} \mid {\mathfrak A}_t \right) \cr
& + \beta{\mathbb E} 
\left(
N^*_{t+1}Q^*_{t + 1} \psi_{x}^g ( D_t, X_t, W_{t+1} )'  \widetilde{MG}_{t+1}\mid {\mathfrak A}_t \right).  
\end{align*} 
$$

$$
\begin{align*}
\widetilde{MG}_t = \kappa(D_t, X_t) + {\mathbb E} \left(
N^*_{t+1}Q^*_{t + 1}  \left(\frac { G_{t+1}}{ G_t}\right)  \widetilde{MG}_{t+1}\mid {\mathfrak A}_t \right)
\end{align*}
$$ 

$$ 
\begin{align*}
& 0 =   \kappa_d( D_t, X_t) 
\cr &+ 
\beta 
{\mathbb E} 
\left(
N^*_{t+1}Q^*_{t + 1} \left(\frac {G_{t+1}}{G_t}\right) \psi_{d'}^x ( D_t, X_t, W_{t+1} )' {\widetilde {MX}} _{t+1} \mid {\mathfrak A}_t \right) \cr
&+  \beta 
{\mathbb E} 
\left(
N^*_{t+1}Q^*_{t + 1}  \psi_{d'}^g ( D_t, X_t, W_{t+1} ) {\widetilde {MG}}_{t+1}  \mid {\mathfrak A}_t \right) .
\end{align*}
$$ 

+++ {"id": "966a389a"}

# 4 Algorithm 

## 4.1 Approach 1
The algorithm initializes $\mu^0, \Upsilon^2_0, \Upsilon^2_1,$ and $\Upsilon^2_2$ at zero.  It uses the first and second order approximations of Schmidt-Grohe and Uribe (2004) and Lombardo and Uhlig (2018) with some modifications.  


### 4.1.1: Step 1
Take as input $\mu^0, \Upsilon^2_0, \Upsilon^2_1,$ and $\Upsilon^2_2$.  Transform the first-order and second-order approximations for the state evolution by using the probability induced by $N_{t+1}^0.$  Under this distribution, $W_{t+1}$ has a conditional distribution with a mean given by $\mu^0$ and covariance matrix $I$.

### 4.1.2: Step 2
When solving for the first-order approximation, an additional term needs to added to the forward equation system to accomodate the recursive utility adjustment: 

$$
\bar{H}^1 \stackrel{\text { def }}{=} \frac{(\rho-1)}{2\left(1-\gamma_o\right)}\left|\mu^o\right|^2 H_{t+1}^0 
$$

where

$$
H_{t+1}^0 =  \psi_1(X_t^0, J_t^0,  X_{t+1}^0, J_{t+1}^0, 0,0)
$$

and is time invariant.

### 4.1.3: Step 3
When solving for the second-order approximation, an additional term needs to be added to the forward equation system to accomodate the recursive utility adjustment:


$$
\begin{aligned} 
{\bar H}^2_t  \stackrel{\text { def }}{=} &\left(1-\gamma_o\right) \Theta_2^1\left(\Upsilon_1^2 X_t^1+\Upsilon_0^2\right) \\
& +(\rho-1) \mu^o \cdot\left(\Upsilon_1^2 X_t^1+\Upsilon_0^2\right) H_{t+1}^0 \\
& + 2 \frac{(\rho-1)}{\left(1-\gamma_o\right)}\left[\Theta_2^1 \mu^0+\frac{1}{2} \mu^0 \cdot \mu^0\left(\Theta_0^1+\Theta_1^1 X_t^1\right)\right] \\
& + \left(\frac{1-\rho}{1-\gamma_o}\right){ }^2\left(\left|\mu^0\right|^2+\frac{1}{4}\left|\mu^0\right|^4\right)
\end{aligned} 
$$

where $\Theta_0^1, \Theta_1^1,$ and $\Theta_2^1$ come from expressing the first-order approximation, $H_{t+1}^1,$ of 

$$
H_{t+1} \stackrel{\text { def }}{=}  \psi_1(X_t, J_t, X_{t+1}, J_{t+1})
$$

as

$$
H_{t+1}^1=\Theta_0^1+\Theta_1^1 X_t^1+\Theta_2^1\left(W_{t+1}-\mu^0\right)
$$

### 4.1.4: Step 4

Compute an upated $\mu^0, \Upsilon^2_0, \Upsilon^2_1,$ and $\Upsilon^2_2$ and return to Step 1.

The algorithm continues until it reaches a numerical convergence.

## 4.2 Approach 2

In this aproach, we use an alternative approximation of $N_{t+1}^*$ given by

$$
\log {\widetilde N}_{t+1} =  \left[(1 - \gamma_o)\left({\widehat V}_{t+1}^1 + {\frac 1 2} {\widehat V}_{t+1}^2 \right) \right] - \log {\mathbb E} \left( \exp\left[(1 - \gamma_o)\left({\widehat V}_{t+1}^1 + {\frac 1 2} {\widehat V}_{t+1}^2 \right) \right] \mid {\mathfrak A}_t \right)
$$

The induced distribution of $W_{t+1}$ is normal, with mean $\tilde{\mu}_t$ and variance $\widetilde{\Sigma}$. This distribution is used in both the first-order and second-order contributions to the approximation. The conditional mean for $\tilde{\mu}_t$ is affine in $X_t^1$. 
<br><br>

+++ {"id": "5d5c1350"}

<!-- # 4 Adjustment cost example


## 4.1 Production
Consider a model recursive utility preferences and an AK technology with adjustment costs. 


$$
{C_t} + {I_t}  = {\mathbf a}K_t\tag{12}
$$
$$
K_{t+1}  = K_t \left[1 + \theta_2 \left({\frac {I_t} {K_t}}\right) \right]^{\theta_1} \exp( - \alpha_k + Z_t - {\frac 1 2} \mid \sigma_k \mid^2  + \sigma_k\cdot W_{t+1})\tag{13}
$$

$$
Z_{t+1} = {\mathbb A} Z_t + \mathbb{B} W_{t+1} \tag{14}
$$
where $Z_{t+1}$ is a scalar. $W_{t+1}$ is a shock vector containing 2 entries.  For the purposes of computation, we divide the first two equations by $K_t$.

## 4.2 FOC on Investment
The consumer's stochastic discount factor is

$$
\frac{S_{t+1}}{S_t} = \beta \left(\frac{V_{t+1}}{R_t}\right)^{1-\gamma} \left(\frac{V_{t+1}}{R_t}\right)^{\rho-1} \left(\frac{C_{t+1}}{C_t}\right)^{-\rho}
$$

The first-order condition for investment is:

$$
MC_t \mathbb{E}\left[\left(\frac{S_{t+1}}{S_t}\right) \frac {MK_{t+1}}{MC_{t+1}}  \left(\theta_1 \theta_2\left[1+\theta_2\left(\frac{I_t}{K_t}\right)\right]^{\theta_1-1} \exp(B_{t+1} - B_t) \right) \mid \mathfrak{A}_t\right] - MC_t  = 0.  
$$

where $MK_{t+1}$ is the co-state associated with the capital evolution and $MC_t$ is the multiplier associated with the production relation. Dividing this first-order condition by $MC_t$ gives:

$$
\mathbb{E}\left[\left(\frac{S_{t+1}}{S_t}\right) \frac {MK_{t+1}}{MC_{t+1}}  \left(\theta_1 \theta_2\left[1+\theta_2\left(\frac{I_t}{K_t}\right)\right]^{\theta_1-1} \exp(B_{t+1} - B_t)\right) \mid \mathfrak{A}_t\right] - 1  = 0,\tag{15}
$$

which is the equation we use in computation.  

The co-state equation satisfies:

$$
\begin{align*}
& MC_t \mathbb{E}\left[\left(\frac{S_{t+1}}{S_t}\right) \left(\frac {MK_{t+1}}{MC_{t+1}} \right) \left[\left( \frac{K_{t+1}}{K_t} \right) - \left(\theta_1 \theta_2\left[1+\theta_2\left(\frac{I_t}{K_t}\right)\right]^{\theta_1-1} \exp(B_{t+1} - B_t)\right)\left(\frac {I_t} {K_t} \right)  \right] \mid {\mathfrak A}_t \right]\cr &  - MK_t + {\mathbf a} MC_t = 0
\end{align*}
$$

Dividing this equation by $MC_t$ gives:

$$
 \mathbb{E}\left[\left(\frac{S_{t+1}}{S_t}\right) \left(\frac {MK_{t+1}}{MC_{t+1}} \right) \left[\left( \frac{K_{t+1}}{K_t} \right) - \left(\theta_1 \theta_2\left[1+\theta_2\left(\frac{I_t}{K_t}\right)\right]^{\theta_1-1}  \exp(B_{t+1} - B_t)\right)\left(\frac {I_t} {K_t} \right)  \right] \mid {\mathfrak A}_t \right] - \frac {MK_t}{MC_t} + {\mathbf a}  = 0,
\tag{16}
$$

which is the equation we use in computation.  

## 4.3 State and jump variables

This formulation has two state variables: $\log K_{t} - \log K_{t-1},$ and $Z_{t}$; and has three jump variables:
$\log C_t - \log K_t,\log I_t - \log K_t$ and $\log MK_t - \log MC_t$.   
<br><br> -->

+++

<!-- # 6 Shock elasticities

This section illustrates how to use the results from Expansion Suite to compute shock elasticities of a stochastic process, for example, consumption growth. -->

+++

<!-- ## 6.1 Compute elasticities

To compute exposure elasticity, we may use the `exposure_elasticity` function. The function takes 6 inputs:

`log_M_growth` : *LinQuadVar*

- Consumption growth process, can be loaded from the `ModelSolution` attribute `gc_tp1`.

`X1_tp1` : *LinQuadVar*

- First order expansion of the state evolution processes, can be loaded from the `ModelSolution` attribute `X1_tp1`.

`X2_tp1` : *LinQuadVar*

- Second order expansion of the state evolution processes, can be loaded from the `ModelSolution` attribute `X2_tp1`.

`T` : *int*
 
- Time horizon (in month).

`shock` : *int*

- The shock component we are interested in, `0 <= shock <= n_W-1`.

`percentile` : *float*

- The percentile of the elasticity we want to compute, `0 <= percentile <= 1`. For example, set `percentile = 0.5` to compute the median.

To compute price elasticity, we may use the `price_elasticity` function. The function takes 7 inputs. In addition to the 6 described above, we also need

`log_S_growth` : *LinQuadVar*

- Approximation of the log stochastic discount factor, which can be computed using attributes of `ModelSolution`. -->

+++

<!-- 
Recall that

$$
\begin{aligned}
\frac{S_{t+1}}{S_t} & =  \beta N_{t+1}^*\left(\frac{C_{t+1}}{C_t}\right)^{-\rho}  \exp \left[(\rho-1)\left(\widehat{V}_{t+1}-\widehat{R}_t\right)\right] \\
& = N_{t+1}^* \exp \left(\widehat{S}_{t+1}-\widehat{S}_t\right)
\end{aligned}
$$

where ${\log N^*_{t+1}}$ is from Section 4.3 in the notes and is approximated by `log_N_tilde` in the model solution:

$$
\log {\widetilde N}_{t+1} =  \left[(1 - \gamma_o)\left({\widehat V}_{t+1}^1 + {\frac 1 2} {\widehat V}_{t+1}^2 \right) \right] - \log {\mathbb E} \left( \exp\left[(1 - \gamma_o)\left({\widehat V}_{t+1}^1 + {\frac 1 2} {\widehat V}_{t+1}^2 \right) \right] \mid {\mathfrak A}_t \right)
$$

and

$$
\widehat{S}_{t+1}-\widehat{S}_t \stackrel{\text { def }}{=}\log\beta-\rho\left(\widehat{C}_{t+1}-\widehat{C}_t\right)+(\rho-1)\left(\widehat{V}_{t+1}-\widehat{R}_{t}\right).
$$

We follow approach 3 in the note (Section 4.3) and approximate $\widehat{S}_{t+1}-\widehat{S}_t$ as

$$
\widehat{S}_{t+1}-\widehat{S}_t \approx\left[\widehat{S}_{t+1}^0-\widehat{S}_t^0\right]+\left[\widehat{S}_{t+1}^1-\widehat{S}_t^1\right]+\frac{1}{2}\left[\widehat{S}_{t+1}^2-\widehat{S}_t^2\right]
$$

where

$$
\begin{aligned}
\widehat{S}_{t+1}^0-\widehat{S}_t^0 
    &= \log \beta-\rho \eta_c^0 \\
\widehat{S}_{t+1}^1-\widehat{S}_t^1
    &=-\rho\left(\widehat{C}_{t+1}^1- \widehat{C}_t^1\right)+(\rho-1)\left(\widehat{V}_{t+1}^1-\widehat{R}_t^1\right) \\
\widehat{S}_{t+1}^2-\widehat{S}_t^2 
    &=-\rho\left(\widehat{C}_{t+1}^2-\widehat{C}_t^2\right)+(\rho-1)\left(\widehat{V}_{t+1}^2-\widehat{R}_{t}^2\right)
\end{aligned}
$$

These are implemented in the function `calc_SDF` below. -->

+++

<!-- Compute the 25-th, 50-th, and 75-th percentile of the exposure and price elasticities for 30 years. -->

+++

<!-- ## 6.2 Plot elasticities

In this section, we plot the exposure and the price elasticities of the consumption growth process computed above. -->

+++

<!-- ### 6.2.1 Exposure elasticity -->

+++

<!-- ### 6.2.1 Price elasticity -->

+++

# Reference

[1] Lombardo, Giovanni and Harald Uhlig. 2018. A Theory of Pruning. *International Economic Review* 59 (4):1825–1836.

[2] Schmitt-Grohe, Stephanie and Mart ́ın Uribe. 2004. Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function. *Journal of Economic Dynamics and Control* 28 (4):755–775.

```{code-cell} ipython3

```
