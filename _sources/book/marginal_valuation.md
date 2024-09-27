---
title_to_header: false
---
# Representing Marginal Valuation

**Authors:**  Lars Peter Hansen and Thomas J. Sargent

**Date:** September 2024
$\newcommand{\eqdef}{\stackrel{\text{def}}{=}}$

%Now we can write vectors using our new command: $ f \eqdef x $.





## Introduction


Partial derivatives of value functions measure marginal valuations and appear in first-order conditions of Markov decision problems. They feature prominently in max-min formulations of robust control problems. They can be used to measure losses from suboptimal choices. For example, the social cost of global warming is the important contributor to calculations of the social cost of carbon emissions often inferred from marginal impacts of fossil fuel emissions on climate indicators measured as potentially uncertain damages to economic opportunities in the future.  By importing insights about stochastic nonlinear impulse response functions and from asset pricing methods for valuing uncertain cash flows, this chapter constructs representations of partial derivatives and decompositions of them by i) time horizon and ii) marginal contributions to future utilities.

## Discrete time

We first consider a discrete-time specification.  As in the previous chapter, 
we start with Markov process

```{math}
:label: evolve
X_{t+1}  = \psi(X_t, W_{t+1}), \\ 
Y_{t+1} - Y_t  = \kappa(X_t, W_{t+1}),
```
where there are $n$ components of $X,$  $Y$ is scalar, and $W$ is $k$ dimensional.  Recall the variational processes studied previously:  

```{math}
:label: var_evolve
\Lambda_{t+1}  = \frac {\partial \psi}{\partial x'} (X_t, W_{t+1}) \Lambda_t \\ 
\Delta_{t+1} - \Delta_t  = \frac {\partial \kappa}{\partial x}(X_t, W_{t+1}) \cdot \Lambda_t .
```


We use stochastic impulse responses to provide an "asset pricing" representation of partial derivatives of a value function with respect to one of the components of $X_0$. Consider a value function that satisfies:
```{math}
:label: value_function
\begin{align} 
V(X_t) + Y_t =  & \exp(-\delta) \mathbb{E} \left[ V(X_{t+1}) + Y_{t+1} \mid X_t \right] \\
 &  + [1 - \exp(-\delta)] \left[ U(X_t) + Y_t \right].
 \end{align}
```
This value function need not coincide with a solution to an optimal control problem.  It could just be the evaluation of some non-optimal decision rule associated with a stochastic equilibrium.  Marginal valuation is still used as part of a local policy analysis.  


Differentiate both sides of this {eq}`value_function` with respect to $X_t$ and $Y_t$ and form dot products with appropriate variational counterparts:

```{math}
:label: value_differentiation
\begin{align}
\frac{\partial V}{\partial x}(X_t) \cdot \Lambda_t + \Delta_t = & \exp(-\delta) \mathbb{E}\left[ \frac{\partial V}{\partial x}(X_{t+1}) \cdot \Lambda_{t+1} + \Delta_{t+1} \mid X_t, \Lambda_t, \Delta_t \right] \\ 
& + [1 - \exp(-\delta)]\left[\frac{\partial U}{\partial x}(X_t) \cdot \Lambda_t + \Delta_t \right] 
\end{align} 
```
View equation {eq}`value_differentiation`  as a stochastic difference  equation and solve it forward for $\frac{\partial V}{\partial x}(X_t) \cdot \Lambda_t + \Delta_t:$ 
```{math}
:label: value_derive
\begin{align*}
& \frac{\partial V}{\partial x}(X_t) \cdot \Lambda_t + \Delta_t \\  & = 
[1 - \exp(-\delta)]\sum_{\tau = 0}^\infty \mathbb{E}\left( \exp(-\tau \delta) \left[\frac{\partial U}{\partial x}(X_{t+\tau}) \cdot \Lambda_{t+\tau} + \Delta_{t + \tau}\right] \mid X_t, \Lambda_t, \Delta_t \right)
\end{align*}
```
Initialize $\Lambda_0 = \mathrm{e}_i$ where $\mathrm{e}_i$ is a coordinate vector with a one in position $i$ and $\Delta_0 = 0.$ We now have represented the partial derivative of the value function as:
```{math}
\frac{\partial V}{\partial x_i}(x) = [1 - \exp(-\delta)]\sum_{t = 0}^\infty \exp(-t \delta) \mathbb{E}\left[ \frac{\partial U}{\partial x}(X_{t}) \cdot \Lambda_{t} + \Delta_t \mid X_0 = x, \Lambda_0 = \mathrm{e}_i, \Delta_0 = 0 \right]
```
This resembles an asset pricing formula in which 
```{math}
\left\{\exp(-\delta t)\begin{bmatrix} \Lambda_t \\ \Delta_t \end{bmatrix} : t \ge 0 \right\}
```
acts as a *vector stochastic discount factor process* and the marginal contribution
```{math}
\left\{ [1 - \exp(-\delta)]\begin{bmatrix}
\frac{\partial U}{\partial x}(X_{t}) \\ \Delta_t \end{bmatrix} : t \ge 0\right\},
```
acts as a  *vector stochastic cash flow process*. The stochastic impulse response tells the marginal state vector response at date $t$ to changing the $i^{th}$ state vector component at date zero, while the vector of marginal cash flows at date $t$ measures the impact on utility of a marginal change in the date $t$ state vector.


````{prf:remark}
Sometimes it is convenient to apply summation by parts:
```{math}
\begin{align*} 
& [1 - \exp(-\delta)]\sum_{\tau = 0}^\infty \mathbb{E}\left( \exp(-\tau \delta) \Delta_{t + \tau} \mid X_t, \Lambda_t, \Delta_t \right) \\ 
& = \sum_{\tau = 1}^\infty \mathbb{E}\left[ \exp(-\tau \delta) \left(\Delta_{t + \tau} - \Delta_{t+\tau -1} \right) \mid X_t, \Lambda_t, \Delta_t \right] + \Delta_t \\ 
& = \sum_{\tau = 1}^\infty \mathbb{E}\left[ \exp(-\tau \delta) \frac{\partial \kappa}{\partial x}(X_{t+\tau-1}, W_{t+\tau})\cdot \Lambda_{t+ \tau} \mid X_t, \Lambda_t, \right] + \Delta_t .
\end{align*}
```
Substituting into {eq}`value_derive` gives:
```{math}
\frac{\partial V}{\partial x}(X_t) \cdot \Lambda_t + \Delta_t = [1 - \exp(-\delta)]\sum_{\tau = 0}^\infty \mathbb{E}\left[ \exp(-\tau \delta) \frac{\partial U}{\partial x}(X_{t+\tau}) \cdot \Lambda_{t+\tau} \mid X_t, \Lambda_t, \right]
 \\
 + \sum_{\tau = 1}^\infty \mathbb{E}\left[ \exp(-\tau \delta) \frac{\partial \kappa}{\partial x}(X_{t+\tau-1}, W_{t+\tau})\cdot \Lambda_{t+ \tau} \mid X_t, \Lambda_t, \right] + \Delta_t.
```
````


## Continuous time

A continuous-time counterpart allows us to draw a distinction between small shocks (Brownian increments) and large shocks (Poisson jumps).  Formally, 
we consider a continuous-time specification with Brown motion shocks, in other words, diffusion dynamics.   We then allow for jumps by treating them  as terminal conditions where we impose continuation values conditioned on a jump taking place.   The possibility of a jump will contribute to the value function computation.  After developing the basic approach we extend the analysis to include robustness in the valuation.  



### Diffusion dynamics
As a part of a more general derivation, we begin with state dynamics modeled as a Markov diffusion:

```{math}
\begin{align*}
dX_t & = \mu(X_t) dt + \sigma(X_t) dW_t \\
dY_t & = \nu(X_t) dt + \varsigma(X_t) \cdot dW_t. 
\end{align*}
```
As for discrete time, these dynamics might or might not be outcomes of an optimization problem.

Using the variational process construction in the previous chapter, recall that 
```{math}
d\Lambda_{t}^i = \left(\Lambda_t\right)'\frac{\partial \mu_i}{\partial x}(X_t) dt + \left({\Lambda_t}\right)'\frac{\partial \sigma_i}{\partial x}(X_t) dW_t.
```

With the appropriate stacking, the drift for the composite process $(X,\Lambda)$ is:

```{math}
:label: variational_mean
\mu^a(x,\lambda) \overset{\text{def}}{=} \begin{bmatrix} \mu(x) \\ \lambda'{\frac {\partial \mu_i} {\partial x} }(x) \\
... \\ \lambda'{\frac {\partial \mu_n} {\partial x} }(x) \end{bmatrix},
```
and the composite matrix coefficient on $dW_t$ is given by

```{math}
:label: variational_diffusion
\sigma^a(x,\lambda) \overset{\text{def}}{=} \begin{bmatrix} \sigma(x) \\  \lambda'\frac {\partial \sigma_1 }{\partial x}(x)\\ ... \\
\lambda' \frac {\partial \sigma_n }{\partial x}(x) \end{bmatrix}.
```

Similarly,  $\Delta$ is  the scalar variational process associated with $Y$ with evolution 

```{math}
d \Delta_t = \Lambda_t \cdot \frac {\partial \nu}{\partial x} (X_t)dt + {\Lambda_t}' \frac {\partial \varsigma}{\partial x'} dW_t 
```

### An initial representation of a partial derivative

Consider the evaluation of discounted utility where the instantaneous contribution is $U(x)$ where $x$ is the realization of a state vector $X_t$. The function $U$ satisfies a Feynman-Kac (FK) equation: 

```{math}
:label: FKequation
\begin{align}
0 = & \delta \left[U(x) + y\right]  - \delta \left[V(x) + y \right] 
\mu(x) \cdot \frac {\partial V}{\partial x}(x) + \nu(x) 
\cr &+ {\frac 1 2 }{\rm trace}\left[\sigma(x)' \frac {\partial^2 V}{\partial x \partial x'}(x) \sigma(x) \right].
\end{align}
```
As in the discrete-time example, we want to represent

```{math}
V_{x_i}(x) = {\frac {\partial V}{\partial x_i}}(x)  
```
as an expected discounted value of a marginal impulse responses of future $X_t$ to a marginal change of the $i^{th}$ coordinate of $x.$ 

By differentiating Feynman-Kac equation {eq}`FKequation` with respect to each coordinate, we obtain a vector of equations one for each state variable. We then form the dot product of this vector system with respect to $m$ to obtain a scalar equation system that is of particular interest. The resulting equation is a Feynman-Kac equation for the scalar function:

```{math}
\lambda  \cdot \frac {\partial V}{\partial x}
```
as established in the Appendix. Given that the equation to be solved involves both $\lambda$ and $x$, this equation uses the diffusion dynamics for the joint process $(X,\Lambda)$.

The solution to this Feynman-Kac equation takes the form of a discounted expected value:

```{math}
:label: value_function_partial
\frac {\partial V}{\partial x}(X_0) \cdot \Lambda_0 + \Delta_0 =  \delta  \int_0^\infty  \exp( - \delta t ) {\mathbb E} \left[
\frac {\partial U}{\partial x} (X_{t}) \cdot \Lambda_{t} + \Delta_t  \mid X_0, \Lambda_0, \Delta_0 \right] dt.   
```
By initializing the state vector $\Lambda_0$ to be a coordinate vector of zeros in all entries but entry $i$ and $\Lambda_0 = 0$, we obtain the formula we want, which gives the partial derivative as a discounted present value using $\delta$ as the discount rate. The contribution, $\Lambda_{t},$ is the marginal response of the date $t$ state vector to marginal change in the $i^{th}$ component of the state vector at date zero. The marginal change in the date $t$ state vector induces marginal reward at date $t$:

```{math}
\delta \frac {\partial U}{\partial x} (X_{t})\cdot \Lambda_{t}  + \Delta_t
```
which provides us with a useful interpretation as an asset price. The process $\Lambda$ gives a vector counterpart to a *stochastic discount factor* process and $\delta \frac {\partial U}{\partial x} (X_{t}) + \Delta_t$ gives the counterpart to a *cash flow* to be valued.

One application of representation {eq}`value_function_partial` computes the discounted impulse response:

```{math}
\delta   \exp( - \delta t ) {\mathbb E} \left[
\frac {\partial U}{\partial x_j} (X_{t}) \cdot \Lambda_{t}  \mid X_0, \Lambda_0, \Delta_0 \right].  
```
for $t \ge 0$ and for $j=1,2,...,n$ along with 

```{math}
\delta   \exp( - \delta t ) {\mathbb E} \left[
 \Delta_t  \mid X_0, \Lambda_0, \Delta_0 \right] 
 ```
for $t \ge 0$ as an intertemporal,  additive decomposition of the marginal valuation of one of the state variables as determined by an initialization of $\Lambda_0.$


````{prf:remark}
Representations similar to {eq}`variational_diffusion` appear in the sensitivity analyses of options prices. See {cite}`Fournieetal:1999`.  
````

### Robustness

We next consider a general class of drift distortions that can help us study model misspecification concerns. We initially explore the consequences of exogenously-specified drift distortion. After that, we show how such a distortion can emerge endogenously as a decision-maker's response to concerns about model misspecifications.

For diffusions, we entertain distortions to the Brownian increment.  Instead $W$ being a multivariate Brownian motion, we allow it to have a drift $H$ under a change in the probability distribution.  We index the alternative probability specifications with their corresponding drift processes $H$.  Locally, 
```{math}
dW_t = H_t dt + dW^H_t
```
where $W^H$ is a Brownian motion under the $H$ probability.  Given that both the distribution parameterize by $H$ and the baseline distribution for the increment are normals with an identity matrix as the local covariance matrix, the local measure of relative entropy is given by the quadratic term: 
```{math}
{\frac 1 2} H_t \cdot H_t .
```
See {cite}`James:1992`, {cite}`AndersonHansenSargent:2003`,  and {cite}`hstw:2006` for further discussions.   

Initially, we introduce an exogenously specified drift distortion process $H$ into the diffusion dynamics:
```{math}
\begin{align*}
&d X_t = \mu(X_t)dt + \sigma(X_t) H \left( \bar{X}_t \right) dt + \sigma(X_t) dW_t^H \\
&d \overline{X}_t = \bar{\mu}\left( \overline{X}_t \right) dt +  \bar{\sigma}\left(\overline{X}_t \right) dW_t^H.
\end{align*}
```


By imitating our earlier analysis, we can associate with this joint system $(X, \overline{X})$ a composite variational process $(\Lambda, \overline{\Lambda}).$ To study endogenous state variable sensitivity, we are especially interested in the $M$ component, not the $\overline{\Lambda}_t$ process. Notice that if we set $\overline{\Lambda}_0 = 0$, then $\overline{\Lambda}_t = 0$ for $t > 0.$

The evolution for the variational process component $\Lambda$ becomes:
```{math}
d\Lambda_{t}^i = \left(\Lambda_t\right)'\left[\frac{\partial \mu_i}{\partial x}(X_t)  + \frac{\partial \sigma_i}{\partial x}(X_t) H\left( \overline{X}_t \right) \right]dt 
+
\left(\Lambda_t\right)'\frac{\partial \sigma_i}{\partial x}(X_t) dW_t^H.
```
Importantly, there is no contribution from differentiating $H$ with respect to $x$ since $H$ only depends on the $\bar{X}_t$ process.

````{prf:remark}
While we focus on Markov forms of misspecification, this can be relaxed.  The misspecification that will  most concern a  decision maker will have a Markov representation. That helps explain why we  make a Markov assumption here.  
````

### Value function derivatives under robustness

We now  let the flow term  be: 
```{math}
\delta \left[ U(x) + y \right] +{\frac \xi 2} \vert H\left( \bar{x} \right) \vert^2.
```
This implies a value function ${\overline V}(x,{\bar x}) + y.$ 

Consider another value function that is sometimes  used to compute a robustness adjustment to valuation.  It coincides with ones sometimes used in robust control problems.   This value function solves the HJB equation:
```{math}
:label: min:hjb
\begin{align}
0 = \min_h & \hspace{.2cm} \delta \left[U(x) + y\right]  - \delta \left[ V(x) + y \right]+ {\frac \xi 2}|h|^2 + \left[\mu(x) +\sigma(x)h \right] V_x(x) \cr
& + \nu(x) + \varsigma(x) \cdot h +
{\frac 1 2} {\rm trace} \left[ \sigma (x)' V_{xx}(x) \sigma(x) \right].
\end{align}
```
The first-order conditions for $h$ in equation {eq}`min:hjb` imply that 
```{math}
\sigma(x)'V_x(x)  + \varsigma(x) + \xi h = 0  
```
The solution, $V,$ of the  HJB equation satisfies the Feynman-Kac equation that emerges after  we substitute the minimizing $h$ into the HJB equation:

```{math}
:label: HJB_robust
\begin{align}
0 =& \hspace{.2cm} \delta U(x) - \delta V(x) + {\frac \xi 2}|h^*(x) |^2 + \left[\mu(x) +\sigma(x)h^*(x)  \right] \cdot  V_x(x) \\ & + \left[\nu(x) + \varsigma(x)\cdot h^*(x) \right]
 + {\frac 1 2} {\rm trace} \left[ \sigma (x)' V_{xx}(x) \sigma(x) \right] ,
\end{align}
```
where 
```{math}
:label: optimal_distortion
h^*(x) = - \frac 1 \xi \left[ \sigma'(x) V_x(x) + \varsigma(x) \right].
```

Consider an exogenously specified drift distortion as in the previous subsection where $H = h^*$ and  stochastic dynamics  $\bar X$ satisfy a consistency requirement:

```{math}
\begin{align}
{\bar \mu }(x)  &= \mu(x) + H(x) \\
{\bar \sigma}(x) &= \sigma(x) 
\end{align}
```
and ${\bar X}_0 = X_0.$    

We  now show that when 
$ H(\bar x) = h^*( \bar x ),$ it is also true that 

```{math}
:label: result_value
\begin{align}
V(x) & = {\overline V}(x,x) \\
V_x(x) & = {\overline V}_x(x,x) 
\end{align}
```

Differentiate equation {eq}`HJB_robust`  with respect to $x$:
 
```{math}
\begin{align}
0 = & - \delta V_x + \delta U_x + V_{xx}\left(\mu +\sigma  h^*  \right) 
+ (\mu_x)'V_x + {\rm{mat}} \left\{ \left(\frac {\partial \sigma_i}  {\partial x} \right) h^*   \right\}'V_x + \frac{\partial \nu}{\partial x} + \frac {\partial \varsigma'} {\partial x} h^*
\\ & + {\frac \partial {\partial x}} \left[{\frac 1 2} {\rm trace} \left( \sigma' V_{xx} \sigma \right) \right]  
\end{align}
```
where $\rm{mat}$ denotes a matrix formed by stacking the column arguments.   (This expression uses  the first-order conditions for h and an "Envelope Theorem" to cancel some terms.)
Take corresponding  derivatives  of ${\overline V}$ with respect to the first-argument and then substitute $x={\bar x}$ to obtain {eq}`result_value` when we set $H(x) = h^*(x).$

Note that it follows from the second equation in {eq}`result_value` that
```{math}
V_{xx}(x) = {\overline V}_{xx}(x,x) + {\overline V}_{x\bar x}(x,x) 
```

````{prf:remark}
We can drop $\frac{\xi}{2} |H(\overline{X}_t)|^2$ from the flow term that we used to construct $\overline{V}$ and still obtain the second equality in {eq}`result_value` involving first-derivatives of value functions.
````

We can now compute representations that we have been seeking by simulating under the endogenously determined worst-case probability specification; we can obtain decompositions of various contributions over time and state vector components based on:

```{math}
:label: value_function_partial_robust
\begin{align*}
&\frac{\partial \overline{V}}{\partial x}(X_0, \overline{X}_0) \cdot \Lambda_0 + \Delta_0 \cr
&= \delta \int_0^\infty \exp(-\delta t) \widetilde{\mathbb{E}} \left[ \frac{\partial U}{\partial x}(X_{t}) \cdot \Lambda_{t} + \Delta_t \mid X_0, \bar{X}_0, \Lambda_0, \Delta_0 \right] dt
\end{align*}
```
where set $\overline{X}_0 = X_0$ and $\Lambda_0$ equal to one of the coordinate vectors and $\Delta_0 = 0.$ The mathematical expectation, $\widetilde{\mathbb{E}}$, is computed under the worst-case stochastic evolution computed by imposing $H(\overline{X}_t) = h^*(X_t).$

````{prf:remark}
While we demonstrated that we can treat a drift distortion as exogenous to the original state dynamics, for some applications we will want to view it as a change in the endogenous dynamics that are reflected {eq}`optimal_distortion`.
````

````{prf:remark}
By construction, along a simulated path, $\overline{X}_t = X_t$, which we impose in our numerical calculations. As a consequence, it suffices to simulate $(X,\Lambda)$.    
````

### Allowing IES to differ from unity

Let $\rho$ be the inverse of the intertemporal elasticity of substitution for a recursive utility specification. The utility recursion is now:

```{math}
\left(\frac{\delta}{1-\rho}\right)\left(\exp\left[(1-\rho)\left[U(X_t)+Y_t\right]\right)\exp\left[(\rho-1)\left[\overline{V}(X_t,\bar{X}_t)+Y_t\right]\right]-1\right) + \overline{\mu}_{v,t} = 0
```

where $\overline{\mu}_{v,t}$ is the local mean of $\overline{V}(X,\bar{X}) + Y$ with the robust adjustment discussed above. Compute:

```{math}
\begin{align*}
&\frac{\partial}{\partial x}\left(\frac{\delta}{1-\rho}\right)\left(\exp\left[(1-\rho)\left[U(x)+y\right]\right)\exp\left[(\rho-1)\left[\overline{V}(x,\overline{x})+y\right]\right]-1\right) \\
&= \delta\exp\left[(1-\rho)\left[U(x)-\overline{V}(x,\overline{x})\right]\right]\left[\frac{\partial U}{\partial x}(x)-\frac{\partial\overline{V}}{\partial x}(x,\overline{x})\right]
\end{align*}
```

With this computation, we modify the previous formulas by replacing the subjective discount factor, $\exp(-t\delta),$ with

```{math}
D_t \eqdef \exp\left(-\int_0^t\delta\exp\left[(1-\rho)\left[U(X_\tau)-\overline{V}(X_\tau,\overline{X}_\tau)d\tau\right]\right]\right).
```
Thus the instantaneous discount rate is now state dependent and depends on the both how the current utility compares to the continuation value and on whether $\rho$ is greater or less than one.  When the current utility exceeds the continuation value, the discount rate is scaled up when $\rho$ exceeds one and is scaled down when $\rho$ is less than one. 
The instantaneous flow term, $\delta\frac{\partial U}{\partial x}(X_t),$ with

```{math}
  \delta\exp\left[(1-\rho)\left[U(X_\tau)-\overline{V}(X_\tau,\overline{X}_\tau)\right]\right]\frac{\partial U}{\partial x}(X_t).
```
Combining these contributions gives: 

```{math}
\begin{align*} 
& \frac {\partial {\overline V}}{\partial x}\left( X_0, {\overline X}_0 \right) + \Delta_0 \cr 
& =   \delta \int_0^\infty \widetilde{\mathbb E} \left[ D_t 
\exp\left[(1-\rho)\left[U(X_t)-{\overline V}\left( X_t, {\overline X}_t \right) \right]\right]\frac{\partial U}{\partial x}(X_t) \mid X_0 \Lambda_0, \Delta_0 \right] dt.
\end{align*} 
```






````{prf:remark}
When we conduct simulations, we can impose that $\overline{X}_t = X_t$ along with {eq}`result_value`, implying that
```{math}
\overline{V}\left(X_\tau,\overline{X}_\tau\right) = V(X_\tau).
```

````




### Jumps

We study a pre-jump functional equation in which jump serves as a continuation value. We allow multiple types of jumps, each with its own state-dependent intensity. We denote the intensity of a jump of type $\ell$ by $\mathcal{J}^\ell(x)$; a corresponding continuation value after a jump of type $\ell$ has occurred is $V^\ell(x)+y$. In applications, we'll compute post-jump continuation value $V^\ell$, as components of a complete model solution.  To simplify the notation, we impose that $\rho - 1,$ but it is straightforward to incorporate the $\rho \ne 1$ extension we discussed in the previous subsection.  

As in {cite}`AndersonHansenSargent:2003`, an HJB equation that adds concerns about robustness to misspecifications of jump intensities as well as diffusion dynamics is:

```{math}
\begin{align*}
0 = \min_{h, \, g^\ell \text{ for } \ell=1,...,L} & - \delta V + \delta U + {\frac{\xi}{2}}|h|^2 +\left[\mu +\sigma h\right]\cdot \frac {\partial V}{\partial x}  + \nu + \varsigma \cdot h\\
& + {\frac{1}{2}}{\rm trace}\left[\sigma'\frac {\partial^2 V }{\partial x \partial x'}\sigma\right] \\
& + \sum_{\ell=1}^L g^\ell\mathcal{J}^\ell(x)\left[V^\ell - V)\right] \\
& + \xi\sum_{\ell=1}^L\mathcal{J}^\ell(x)\left[1 - g^\ell + g^{\ell}\log g^\ell\right],
\end{align*}
```
where $g^\ell \ge 0$ alters the intensity of type $\ell.$ The term 
```{math}
\sum_{\ell=1}^L\mathcal{J}^\ell(x)\left[1 - g^\ell + g^{\ell}\log g^\ell\right]
```
measures the relative entropy of the alternative jump intensity specifications.  

First-order conditions for the $g^\ell$s are

```{math}
\left[V^\ell(x) - V(x)\right] + \xi\log g^{\ell} = 0.
```

First-order conditions for $h$ remain the same as before.

We again construct a Feynman-Kac equation by substituting in $h^*(x)$ and $g^{\ell*}(x)$. Applying an Envelope Theorem to first-order conditions for minimization tells us that neither $h^*(x)$ nor $g^{\ell*}(x)$ contribute to the derivatives of the value function. This leads us to:

```{math}
:label: HJB_robustjump
\begin{align*}
0 = & -\delta \frac {\partial V }{\partial x} + \delta  \frac {\partial U }{\partial x} +  \frac {\partial^2 V }{\partial x \partial x'}\left(\mu +\sigma h^*\right)\\
& + \left( \frac {\partial \mu'}{\partial x} \right) \frac {\partial V }{\partial x} + {\rm{mat}}\left\{\left(\frac{\partial \sigma_i}{\partial x}\right)h^*\right\}' \frac {\partial V }{\partial x}\\
& +\frac{\partial}{\partial x}\left[\frac{1}{2}{\rm trace}\left(\sigma'  \frac {\partial^2 V }{\partial x \partial x'} \sigma\right)\right] \\
& +\sum_{\ell=1}^L\frac {\partial \mathcal{J}^{\ell}}{\partial x}  g^{\ell*}\left[V^\ell - V\right] \\
& +\sum_{\ell=1}^L\mathcal{J}^{\ell}g^{\ell*}\left[\frac {\partial V^\ell }{\partial x} - \frac {\partial V }{\partial x}\right] \\
& +\xi\sum_{\ell=1}^L\frac {\partial \mathcal{J}^{\ell}}{\partial x} \left[1 - g^{\ell*} + g^{\ell*}\log g^{\ell*}\right].
\end{align*}
```

Applying our $(X,\overline{X})$ analysis tells us that the date $t$ intensity distortion constructed with the minimizing $g^{\ell*}$ depends on the exogenous state vector $\overline{X}_t$ rather than on $X_t$. An Envelope Theorem is again at play here.

It is revealing  to rewrite equation {eq}`HJB_robustjump` as:

```{math}
\begin{align*}
0 = & -\left(\delta + \sum_{\ell=1}^L\mathcal{J}^{\ell}g^{\ell*}\right)\frac {\partial V }{\partial x} + \delta \frac {\partial U }{\partial x} \\
& + \frac {\partial^2 V }{\partial x \partial x'}\left(\mu +\sigma h^*\right) + \left( \frac {\partial \mu'}{\partial x}\right)'\frac {\partial V }{\partial x} + {\rm{mat}}\left\{\left(\frac{\partial \sigma_i}{\partial x}\right)h^*\right\}'\frac {\partial V }{\partial x} \\
& + \frac{\partial}{\partial x}\left[\frac{1}{2}{\rm trace}\left(\sigma'\frac {\partial^2 V }{\partial x \partial x'}\sigma\right)\right] \\
& + \sum_{\ell=1}^L\frac {\partial \mathcal{J}^{\ell}}{\partial x}  g^{\ell*}\left(V^\ell - V\right) \\
& + \sum_{\ell=1}^L\mathcal{J}^{\ell}g^{\ell*} \frac {\partial V^\ell}{\partial x}  \\
& + \xi\sum_{\ell=1}^L \frac {\partial \mathcal{J}^{\ell}}{\partial x}\left[1 - g^{\ell*} + g^{\ell*}\log g^{\ell*}\right].
\end{align*}
```

Notice how distorted intensities act like endogenous discount factors in this equation. The last three terms add a flow term to pertinent Feynman-Kac equations via dot products with respect to $m$. It is significant that these terms do not include derivatives of $g^{\ell*}$ with respect to $x$.

For simulating our asset pricing representation of the partial derivatives of the value function, the discounting term becomes state dependent in order to adjust for the jump probabilitites:

```{math}
D_t \eqdef \exp\left( - \int_0^t\left[\delta +  \sum_{\ell=1}^L\mathcal{J}^{\ell}(X_u)g^{\ell*}(X_u)\right]du\right),
```

In addition, four flow terms are discounted:

```{math}
:label: robustjump_contributions
\begin{align} 
& \delta \Lambda_t \cdot U_x(X_t) & \text{i)}\\
& + \Lambda_t \cdot \sum_{\ell=1}^L\mathcal{J}^{\ell}_x(X_t)g^{\ell*}(X_t)\left[V^\ell(X_t) - V(X_t)\right] &  \text{ii)}\\
& + \Lambda_t \cdot \sum_{\ell=1}^L\mathcal{J}^{\ell}(X_t)g^{\ell*}(X_t) {V^\ell_x}(X_t) & \text{iii)} \\
& + \xi \Lambda_t \cdot\sum_{\ell=1}^L  \mathcal{J}^\ell_x (X_t)\left[1 - g^{\ell*}(X_t) + g^{\ell*}(X_t)\log g^{\ell*}(X_t)\right]. & \text{iv)}
\end{align}
```
Its revealing to think of right side as providing four different sources of the marginal values.  The contributions of $V^{\ell*}$ and $V^\ell_x$ are to be expected because they help to quantify the consequences of potential jumps.  The contribution of $V$ is present because changes in a future state at a given date alters the jump possibility at subsequent dates.  Finally, the fourth term is present because changing a future state changes the penalty for exploring alternative probability distributions.  Increasing this penalty limits the extent of the uncertainty adjustment due to robustness concerns.   

Notice that terms ii) and iv) of formula {eq}`robustjump_contributions` include derivatives of the jump intensity with respect to the state of interest.  In some examples, the jump intensities are constant or depend only on an exogenous state.  I this cases both of the term drop out and only the first and third remain.  


Simulation-based methods can be used to compute these value contributions.  They  should be conducted under implied worst-case diffusion dynamics. With multiple jump components, we decompose contributions to the marginal utility by jump types $\ell$.

### Stochastic growth example

We revisit a model of {cite}`CagettiHansenSargentWilliams:2002` with a finite number of uncertain growth rates. [^hamilton]  

```{math}
:label: capital_evolution
\begin{equation}
dK_t = \left[ \left(Y_t\right)^{1-\alpha}\left(K_t\right)^\alpha - C_t - \eta K_t \right]dt  + K_t \sigma_K \cdot dW_t
\end{equation}
```
where $K_t$ is the date $t$ stock of capital, $C_t$ is consumption, $Y_t$ is an exogenous labor-augmenting technology, labor is in fixed supply, $\alpha$ the relative share of capital in the production of output,  $\eta$ is a depreciation parameter and $\sigma_k \cdot dW_t$ governs the locally stochastic return to capita.  The technology process evolves as a jump process, where
```{math}
:label: productivity_evolution
\begin{equation}
dY_t = Y_t \left(Z_t \cdot \kappa\right)dt   + Y_t \sigma_Y \cdot dW_t, 
\end{equation}
```
and $Z$ evolves as a finite-state Markov chain with a constant intensity matrix $\Lambda$, where a realization of $Z_t$ is on the coordinate vectors in ${\mathbb R}^L,$   We denote the intensity of going from state $j$ to state $\ell$ by  ${\mathcal J}^{j,\ell},$ which is also    the $(j,\ell)$ element of $\Lambda$ is ${\mathcal J}^{j,\ell}$ for $\ell \ne j.$ The diagonal entry, $j,$  of $\Lambda$ is given by :
```{math} 
- \sum_{\ell = 1, \ell  \ne j}^L {\mathcal J}^{j,\ell} . 
```

[^hamilton]: This incorporates an uncertain specification of growth that was originally advocated by {cite}`Hamilton:1998` in his study of output growth.

Transform the two state evolutions  by taking logarithms:

```{math}
\begin{align*}
d \log  K_t & = \left[ \left(\frac {Y_t}{K_t} \right)^{1-\alpha} - \frac {C_t}{K_t}  - \eta  - {\frac 1 2} \vert \sigma_K \vert^2 \right]dt  + \sigma_K \cdot dW_t \cr
d \log Y_t &  = \left[ \left(Z_t \cdot \kappa\right)dt  - {\frac 1 2} \vert \sigma_Y \vert^2\right] dt + \sigma_Y \cdot dW_t,
\end{align*}
```
where we have included local variance adjustments in accordance to Ito's formula.  

For this model, it is convenient  to use $X \eqdef \log K - \log Y$ as the endogenous state variable along with $\log K$, $Z$ the exogenously specified discrete-state jump process, and a control $U \eqdef  \frac {C} {K}.$ This construction of $X$ allows us to impose a particular structure on the value functions, one that simplifies the computations.   The state evolutions that we use for the determination of value functions are:
```{math}
\begin{align*}
d \log K_t = &  \left( \exp\left[(1-\alpha) X_t \right] - U_t - \eta - {\frac 1 2} \vert \sigma_K \vert^2 \right)dt + \sigma_K \cdot dW_t, \cr
d X_t  =  & \left( \exp\left[(1-\alpha) X_t \right]   - U_t  -  \left(Z_t \cdot \kappa\right) - \eta  - {\frac 1 2} \vert \sigma_K \vert^2 + {\frac 1 2} \vert \sigma_Y \vert^2 \right)dt \cr &  + \sigma_X \cdot dW_t 
\end{align*}
```
where $\sigma_X \eqdef  \sigma_K - \sigma_Y.$

This example is under construction with full and partial information. 



### Climate change example

{cite}`BarnettBrockHansenZhang:2024` use representations {eq}`robustjump_contributions` to decompose their model-based measure of the social cost of climate change and the social value of research and development.  In their analysis, there are two types of Poisson jumps.  One is the discovery of a new technology and the other is recognization of how curved the damage function is for more extreme changes in temperature.  The magnitude of damage curvature is revealed by a jump triggered by a temperature anomaly between 1.5 and 2 degrees celsius. {cite}`BarnettBrockHansenZhang:2024` allow for twenty different damage curves. 
While there are twenty one possible jump types, we group them into damage jumps (one through twenty) and a technology jump (twenty one). {cite}`BarnettBrockHansenZhang:2024`
 display the quantitative importance of a technology jump and a damage jump in contributing to the social value of research and development.  We report analogous findings for the social cost of climate change measured as the negative of the marginal value of temperature.  We take the negative of the marginal value of climate change because warming induces a social cost (a negative benefit).

[Table 1](#table1) reports the the contribution of each of the four contributions given on the right side of 
{eq}`robustjump_contributions` for the   partial derivative of the value function with respect to the temperature state variable. The column ii(dc) includes only the contributions for damage curve jumps and column ii(td) includes the remaining contribution from the technology jump, and similarly for iii(dc) and iii(td).  We see from this table that the second term dominates the calculation.  Moreover, the primary source of this contribution is the damage curve uncertainty.  The fourth term partially offsets the cost because of the uncertainty adjustment.  We report this decomposition for three different values of $\xi.$  Decreasing $\xi$ increases the misspecification aversion.  Setting $\xi = \infty$ eliminates robustness concerns altogether.  




<div style="font-size: 16px;">

<!-- <table>
    <thead>
        <tr>
            <th>Case</th>
            <th>i</th>
            <th>ii</th>
            <th>iii(dc)</th>
            <th>iii(td)</th>
            <th>iv</th>
            <th>Sum</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>&#958; = 10000</td>
            <td>-0.0017</td>
            <td>-0.0583</td>
            <td>-0.0009</td>
            <td>-0.0031</td>
            <td>0.0000</td>
            <td>-0.0639</td>
        </tr>
        <tr>
            <td>&#958; = 0.150</td>
            <td>-0.0036</td>
            <td>-0.1387</td>
            <td>-0.0004</td>
            <td>-0.0044</td>
            <td>0.0290</td>
            <td>-0.1180</td>
        </tr>
        <tr>
            <td>&#958; = 0.075</td>
            <td>-0.0070</td>
            <td>-0.4026</td>
            <td>0.0004</td>
            <td>-0.0042</td>
            <td>0.1336</td>
            <td>-0.3382</td>
        </tr>
    </tbody>
</table> -->

<p><strong>Table 1:</strong> 

</div>

| $\xi$         | i       | ii(dc)  | ii(td)  | iii(dc) | iii(td) | iv      | sum     |
|---------------|---------|---------|---------|---------|---------|---------|---------|
| $\infty$      | 0.002   | 0.062   | -0.003  | 0.001   | 0.003   | -0.000  | 0.064   |
| $.150$       | 0.004   | 0.143   | -0.005  | 0.000   | 0.004   | -0.029  | 0.118   |
| $.075$       | 0.007   | 0.408   | -0.004  | -0.000  | 0.004   | -0.134  | 0.280   |

**Table 1:** Components to the partial derivative of the value function with respect to the temperature state variable. The (dc) columns include only the contributions from the damage curve realization jumps, and the (td) columns include the contributions from the technology discovery jump. 


{numref}`discounted_value_of_term2` reports the contributions to the social cost of climate change by horizon for the alternative values of the robustness parameter $\xi$.  Lowering $\xi$ increases the cost contributions at all horizons.  For all specifications of $\xi,$ there is a substantial peak at around 25 years.  To explore further the peak response report densities that capture the timing of the first Poisson jump in {numref}`jump_densities`. We see that these densities peak at around a thirty-year horizon.  Apparently the expected cost contributions peak earlier because for their forward-looking nature.   


```{figure} ../images/stochastic_figure3.png
:name: discounted_value_of_term2
:width: 500px
An intertemporal decomposition of the social cost of temperature for alternative values of $\xi$.  We report only the cost contributions from the second term in {eq}`robustjump_contributions`.   
```





```{figure} ../images/stochastic_figure2.png
:name: jump_densities
:width: 500px
Densities for the  first jump for alternative  values of $\xi$.
```




<!-- ### References

```{bibliography}
:filter: docname in docnames
```

### Footnote -->




<!-- `Bowen or Chun Hei':  please comment out what follows but do not delete.  Relabel the file marginal_valuation and add a new chapter just after stochastic responses called: Representing Marginal Valuation.  Thanks.



Let ${\mathfrak F}$ denote the filtration inclusive of the growth state indicator, $Z,$ and let ${\overline  {\mathfrak F}}$ denote the smaller filtration exclusive of the growth state.  In the latter case, the investors solve a filtering problem to construct the best estimate or conditional expectation, ${\overline Z}_t$,  of the $Z_t$ given currently available information.  Notice the conditional expectation vector, ${\overline Z}_t,$ is the vector of conditional probabilities of each of the $L$ states.



The solution to this model is captured by  coupled Hamilton-Jacobi-Bellman equations. Let $V^j(x)+ \log k $ denote the value function for growth state $j$ where we have  including a conjecture about  dependence of  $V^j$ onto  $\log k$.   Here we are using lower case letters to denote hypothetical realizations of the state variables.  The equations to be solved are 

```{math}
\begin{align*}
0 = \max_u \min_{g,h} & \hspace{.2cm}  \delta \log u   - \delta V^j(x)  + \left([1 + V_x^j(x)]\left[ \exp\left[(1-\alpha) x \right] - u - \eta 
- {\frac 1 2} \vert \sigma_K \vert^2 \right] \right) + \sigma_K \cdot h \cr
& +   V_x^j(x)  \left[  -  z_j  \cdot \kappa +  \sigma_X \cdot h
+ {\frac 1 2}\vert\sigma_Y \vert^2  \right]   \cr 
& +{\frac 1 2} \vert \sigma_X \vert^2 V_{xx}^j(x) 
+ \sum_{\ell =1, \ell \ne j}^L  {\mathcal J}^{j,\ell} g^{\ell}  \left[V^\ell(x) - V^j(x) \right] \cr
& + \xi \frac {h'h} 2 + \xi  \sum_{\ell =1, \ell \ne j}^L {\mathcal J}^{j,\ell} \left( 1 - g^j + g^j \log g^j \right) 
\end{align*} 
```
for $V^j, j=1,2,...,L,$ where $z_j$ is a coordinate vector with a one in position $j$.  For this example, terms ii) and iv) of {eq}`robustjump_contributions` are zero because, conditioned the current growth state, the intensities are independent of $\log K$ and $X$.  



Consider next the case in which the discrete growth states are disguised from investors in the market.   Let ${\mathfrak F}$ be the original filtration (information structure), inclusive of the hidden states,  and let $\overline{\mathfrak F}$ be the smaller filtration that excludes the hidden states.  We use a filtering algorithm derived by {cite}`Wonham:1964` and are motived in part by prior research of {cite}`David:1997` and {cite}`Veronesi:2000`.  To accommodate hidden states and robustness in continuous time, we follow In what follows, we draw on the formulation of {cite}`HansenMiao:2022`.  


We consider the evolution of $\log Y$ and continuation values under the two information structures,   
```{math}
\begin{align*}
d Y_t & =  \left[ \left(Z_t \cdot \kappa\right)dt  - {\frac 1 2} \vert \sigma_Y \vert^2\right] dt + \sigma_Y \cdot  d {W}_t \cr 
d Y_t & =  \left[ \left({\overline Z}_t \cdot \kappa\right)dt  - {\frac 1 2} \vert \sigma_Y \vert^2\right] dt +  \sigma_Y \cdot W_t
\end{align*} 
```
where 
```{math}
d{\overline W}_t \eqdef \begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1}   \begin{bmatrix}   0 \cr
\left(Z_t - {\overline Z}_t \right)  \cdot \kappa \end{bmatrix} dt + dW_t 
```
is a Brownian motion under the filtration $\overline{\mathfrak F}$.  Moreover,
```{math}
\begin{align*}
d {\overline Z}_t & =  \Lambda' {\overline Z}_t dt +   \sigma_{{\bar z}} \left( {\overline Z}_t \right ) d {\overline W}_t \cr
=& \Lambda' {\overline Z}_t dt +   \sigma_{{\bar z}} \left( {\overline Z}_t \right )\begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1}   \begin{bmatrix}   0 \cr
\left(Z_t - {\overline Z}_t \right)  \cdot \kappa \end{bmatrix} dt \cr
& +  \sigma_{{\bar z}} \left( {\overline Z}_t \right ) dW_t 
\end{align*} 
 ```
 where
 ```{math}
 \sigma_{{\bar z}} \left( {\overline Z}_t \right ) \eqdef \left[ \text{diag} \left( {\overline  Z}_t \right) -  {\overline Z}_t {{\overline Z}_t}' \right] \begin{bmatrix} 0 & \kappa   \end{bmatrix} \begin{bmatrix} {\sigma_K} &  {\sigma_X}\end{bmatrix}^{-1} .
 ```
 
For pedagogical reasons, we write the updating equation for prediction using observables as;
```{math}
\begin{align*}
& d {\overline Z}_t =  \Lambda' {\overline Z}_t dt \cr & +  \left[ \text{diag} \left( {\overline  Z}_t \right) -  {\overline Z}_t {{\overline Z}_t}' \right] \begin{bmatrix} 0 & \kappa   \end{bmatrix} \left[\begin{bmatrix} {\sigma_K}' \cr {\sigma_X}'\end{bmatrix}
\begin{bmatrix} {\sigma_K} & {\sigma_X}\end{bmatrix}\right]^{-1} \begin{bmatrix} d \log K_t -  {\bar \mu}_t^K \cr
dX_t - {\bar \mu}_t^X \end{bmatrix}, 
\end{align*}
```
which shows  the contribution that the observation vector, net of its conditional mean, contributes to the updated hidden state information.
In this updating equation, we define:
```{math}
\begin{align*}
& {\bar \mu}_t^K \eqdef \left( \exp\left[(1-\alpha) X_t \right] - U_t - \eta - {\frac 1 2} \vert \sigma_K \vert^2 \right)dt 
\cr & {\bar \mu}_t^X \eqdef \left(    {\overline Z}_t \cdot \kappa  - \exp\left[(1-\alpha) X_t \right]   + U_t + \eta  + {\frac 1 2} \vert \sigma_K \vert^2 - {\frac 1 2} \vert \sigma_Y \vert^2 \right)dt, 
\end{align*}
 ``` 

The associated continuation value evolve as[^backwardsde]
```{math}
\begin{align*}
dV_t & =  {\mu}_t^Vdt   +  {\sigma}_t^V \cdot d {W}_t \cr 
dV_t & = {\bar \mu}_t^Vdt+  \sigma_t^V \cdot d {\overline W}_t
\end{align*} 
```
where the first evolution is under full information and the second one is under the hidden state information.  The drift  coefficients are related via:
```{math}
\mu_t^V - {\bar \mu}_t^V   =  {{\sigma}_t^V}' \begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1}   \begin{bmatrix}   0 \cr
\left({\overline Z}_t -  Z_t \right)  \cdot \kappa \end{bmatrix} . 
```

[^backwardsde]:Value functions are forward-looking and solved by backward induction.  Formally, in continuous time they solve backward stochastic differential equations. 

Since the continuation values that interest us come from a value function, $V\left(x, {\bar z}\right) + \log k,$ 
```{math}
\sigma_t^V =    \begin{bmatrix} \sigma_K & \sigma_X &  \sigma_{{\bar z}} \left( {\overline Z}_t \right )' \end{bmatrix}  \begin{bmatrix} 1  \cr V_x \left( X_t, {\overline  Z}_t \right)  \cr V_{\bar z}\left( X_t, {\overline  Z}_t \right)   \end{bmatrix} 
```
Thus 
```{math}
\begin{align*}
\mu_t^V - {\bar \mu}_t^V & = \begin{bmatrix} 1 &  V_x\left(X_t, {\overline Z}_t \right)  & {V_{\bar z}}\left(X_t, {\overline Z}_t\right)  ' \end{bmatrix}  \begin{bmatrix} {\mathbb I} \cr  \sigma_{{\bar z}} \left( {\overline Z}_t \right ) \begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1} \end{bmatrix} 
\begin{bmatrix}   0 \cr
\left(Z_t - {\overline Z}_t \right)  \cdot \kappa \end{bmatrix} \cr &  = \left[
V_x \left( X_t, {\overline Z}_t \right)  + {V_{\bar z}}\left( X_t, {\overline Z}_t \right)' \sigma_{{\bar z}} \left( {\overline Z}_t \right ) \begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1} \begin{bmatrix}   0 \cr 1
\end{bmatrix} \right] \left({\overline Z}_t -  Z_t \right)  \cdot \kappa . 
\end{align*} 
```


Under hidden information, we  replace the robust adjustment for the jump intensity with a robust adjustment to the constructed conditional probabilities of the hidden states along with a full information, robust adjustment to the drift of the full information Brownian motion.  Taken together, these adjustments are determined by solving: 
```{math}
\min_{H, G, {\mathbb E} \left(G \vert  {\overline  {\mathfrak F}}_t \right)=1 } {\mathbb E} \left[ G \left( \mu^V_t - {\bar \mu}^V_t \right)  + \xi_2 G \log G \vert {\overline  {\mathfrak F}}_t \right] + \sigma_t^V \cdot H + \frac {\xi_1} 2 H'H.  
```
In this minimization problem, we introduce two penalty parameters, one for contribution of a drift to the underlying Brownian motion and another for altering the hidden state probabilities.  While the drift distortion could depend on the hidden states, the minimizing solution does not depend on these states.  


The first-order conditions for minimization imply that 
```{math}
\begin{align*}
& G^*_t = \frac { \exp\left[- \frac 1 {\xi_2} \left( \mu^V_t - {\bar \mu}^V_t \right) \right] } 
{{\mathbb E}\left(  \exp\left[- \frac 1 {\xi_2} \left( \mu^V_t - {\bar \mu}^V_t \right) \right]  \mid  {\overline  {\mathfrak F}}_t \right),} \cr
& H^*_t  = - \left(\frac 1 {\xi_2} \right) \sigma_t^V.
\end{align*}
```
The resulting HJB equation with these robust adjustments is 
```{math}
\begin{align*}
0 = \max_u \min_{h, g, g \cdot {\bar z}  = 1  } & \hspace{.2cm}  \delta \log u   - \delta V(x,{\bar z})  \cr & + \left([1 + V_x(x, {\bar z} )]\left[ \exp\left[(1-\alpha) x \right] - u - \eta
- {\frac 1 2} \vert \sigma_K \vert^2 \right] \right) + \sigma_K \cdot h \cr
& +   V_x(x,{\bar z} )  \left[ -  {\bar z} \cdot \kappa  +  \sigma_X \cdot h
+ {\frac 1 2}\vert\sigma_Y \vert^2  \right]   \cr 
& + V_{{\bar z}} (x, {\bar z})' \Lambda' {\bar z} \cr
&  + \text{trace} \left\{\begin{bmatrix} {\sigma_X}' \cr  \sigma_{\bar z}({\bar z})'\end{bmatrix} \begin{bmatrix} V_{xx} (x, {\bar z}) & V_{x{\bar z}'} (x, {\bar z} ) \cr V_{{\bar z} x} (x, {\bar z} ) &   V_{{\bar z} x}(x, {\bar z}) \end{bmatrix} \begin{bmatrix} {\sigma_X} & \sigma_{\bar z}({\bar z}) \end{bmatrix}  \right\}\cr 
& + \sum_{\ell = 1}^L   \left[
V_x \left( x, {\bar z}  \right)  + {V_{\bar z}}\left( x, {\bar z}  \right)' \sigma_{{\bar z}} \left( {\bar z}  \right ) \begin{bmatrix}   {\sigma_K}' \cr {\sigma_X}' \end{bmatrix}^{-1} \begin{bmatrix}   0 \cr 1
\end{bmatrix} \right] \sum_{\ell = 1}^L  g_\ell {\bar z}_\ell \left(  {\bar z}  \cdot \kappa  - \kappa_\ell  \right)  \cr
& + \xi_1 \frac {h'h} 2 + \xi_2\sum_{\ell=1}^L   g_{\ell} {\bar z}_{\ell}  \log g_\ell . 
\end{align*} 
``` ]:#

 -->
