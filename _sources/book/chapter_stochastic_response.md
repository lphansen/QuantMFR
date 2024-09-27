---
title_to_header: false
---
# Stochastic Responses 

**Authors:** Jaroslav Borovicka, Lars Peter Hansen, and Thomas J. Sargent

**Date:** September 2024
$\newcommand{\eqdef}{\stackrel{\text{def}}{=}}$

%Now we can write vectors using our new command: $ f \eqdef x $.





## Introduction


Impulse-response methods have been used by economists since {cite}`frisch33` and in other disciplines. For nonlinear stochastic models, impulse responses are themselves stochastic.  Alternative approaches have been suggested in economics including {cite}`GallantRossiTauchen:1993`,  {cite}`KoopPesaranPotter:1996`, and {cite}`GourierouxJasiak:2005`.  This chapter provides  stochastic responses in both discrete and continuous time for marginal changes in state variables and shocks.  The vast literature on vector autoregressions views these responses as ends in and of themselves.  Identified shocks as exogenous inputs into a dynamical system in effect `cause` movements in the vector time series of interest. Such measurements, while interesting, can have rather indirect connections to hypothetical interventions related to perspective policy changes that are at the heart of structural models.    Builders of dynamic stochastic equilibrium models use the construct of 'structural' in the sense  of {cite}`Marschak:1953`, {cite}`Hurwicz:1966`, and {cite}`lucas`.   They allow for investigating how a dynamical system changes when one portion of it is altered.  For us, these stochastic responses  are central inputs into marginal valuations that we will use for a variety of purposes.   In subsequent chapters, we provide extensive discussions of two types of such applications.     The first type 
provides asset-pricing type representations for endogenous  variables including various forms of capital.  These include policy relevant variables such as the social cost of climate change and the social value of research and development.   The second type generates what we call shock elasticities that help us characterize the building  blocks for exposures to uncertainty and prices of those exposures.

---


[Partial derivatives of value functions measure marginal valuations and appear in first-order conditions of Markov decision problems. They feature prominently in max-min formulations of robust control problems. They can be used to measure losses from suboptimal choices. For example, the social cost of global warming is the important contributor to calculations of the social cost of carbon emissions often inferred from marginal impacts of fossil fuel emissions on climate indicators measured as potentially uncertain damages to economic opportunities in the future.  By importing insights about stochastic nonlinear impulse response functions and from asset pricing methods for valuing uncertain cash flows, this chapter constructs representations of partial derivatives and decompositions of them by i) time horizon and ii) marginal contributions to future utilities.]:#

## Discrete time

We first consider a discrete-time specification.  

### Discrete-time Markov dynamics

We start with Markov process

```{math}
:label: evolve
X_{t+1}  = \psi(X_t, W_{t+1}), \\ 
Y_{t+1} - Y_t  = \kappa(X_t, W_{t+1}),
```
where there are $n$ components of $X,$  $Y$ is scalar, and $W$ is $k$ dimensional 


### Discrete-time  variational dynamics 

Let $\Lambda$ denote the first variational process for $X$, and let $\Delta$ denote the first variational process for $Y$. These variational processes are the ingredients to stochastic impulse responses to small changes in the underlying state variables. We compute them by "differentiating" in a generalized sense that accommodates the underling stochastic structure.  To obtain a recursive representation for $(\Lambda, \Delta),$  we by differentiate {eq}`evolve` and apply the chain rule:

```{math}
:label: var_evolve
\Lambda_{t+1}  = \frac {\partial \psi}{\partial x'} (X_t, W_{t+1}) \Lambda_t \\ 
\Delta_{t+1} - \Delta_t  = \frac {\partial \kappa}{\partial x}(X_t, W_{t+1}) \cdot \Lambda_t .
```
In this calculation, $\Lambda_{t+1}$ and $\Delta_{t+1}$ are stochastic as they inherit the stochastic dependence of $X_{t+1}$ and $Y_{t+1}.$  By differentiating the process at a given calendar date, we are allowing for date $t$ variables to change as a function of date $t$ information.    

To obtain alternative stochastic (local) impulse response functions, we initialize $({\Lambda_0}', \Delta_0)'$ to be one of the coordinate vectors that depends on one of the initial states that we want to perturb. Then $({\Lambda_t}', \Delta_t)'$ is the date $t$ state vector stochastic response to the perturbation of the initial value of the component.    

To perturb $Y_0,$ we can set $\Lambda_0 = 0$ and $\Delta_0=1;$ then $\Lambda_t = 0$ and $\Delta_t = 1$ for all $t$. Alternatively, if we initialize $\Lambda_0$ be a coordinate vector and $\Delta_0 =0,$ then the response $\Delta_t$ will be a stochastic process.  The outcome of the coordinate vector initialization is a stochastic local impulse response to a marginal change in a particular state variable.  

````{prf:remark}
The evolution of the variational processes is nonstochastic if $ \frac {\partial \psi}{\partial x'}$ and $ \frac {\partial \kappa}{\partial x}$ are constant as  is true when $\psi$ and $\kappa$ are affine in $x$. Otherwise, variational processes are stochastic.
````

````{prf:example} 
:label: ex:quadratic
Consider the following quadratic specification:
```{math}
\begin{align*}
X_{t+1}^i & =  {\sf a}_i \cdot X_t+ {\frac 1 2} {X_t}'{\mathbb A}_i X_t + {X_t}'{\mathbb B}_i W_{t+1} + {\sf b}_i \cdot W_{t+1}, \hspace{.3cm}  i=1,...,n\cr 
Y_{t+1} - Y_t & = {\sf d} \cdot X_t + {\frac 1 2}  {X_t}'{\mathbb D} X_t + {X_t}' {\mathbb F}W_{t+1} + {\sf f}  \cdot W_{t+1} .
\end{align*}
```   
where ${\mathbb A}_i$ and ${\mathbb D}$ are symmetric.  A simple calculation shows:
```{math}
\begin{align*}
\Lambda_{t+1}^i & =  {\sf a}_i \cdot \Lambda_t +  {\Lambda_t}'{\mathbb A}_i X_t + {\Lambda_t}'{\mathbb B}_i W_{t+1} , \hspace{.3cm} i=1,...,n\cr 
\Delta_{t+1} - \Delta_t & =  {\sf d} \cdot \Lambda_t  + {\Lambda_t}'{\mathbb D} X_t + {\Lambda_t}' {\mathbb F}W_{t+1}  .
\end{align*}
```
````




## Continuous-time dynamics

We now consider the continuous-time counterpart for Brownian motion shocks.

### Markov diffusion dynamics 
As a part of a more general derivation, we begin with state dynamics modeled as a Markov diffusion:

```{math}
\begin{align*}
dX_t & = \mu(X_t) dt + \sigma(X_t) dW_t \\
dY_t & = \nu(X_t) dt + \varsigma(X_t) \cdot dW_t. 
\end{align*}
```
where $W$ is now a $k$-dimensional standard Brownian motion.  We denote the filtration (family of specifications of conditioning information events) ${\mathfrak F} \eqdef \left\{ {\mathfrak F}_t  : t\ge 0\right\}$ constructed from the Brownian motion and any pertinent date zero information.  

### Variational process

Following {cite}`BorovickaHansenScheinkman:2014`, we construct marginal impulse response functions using what are called variational processes.  We build the dynamics for what is called the  first valuation  process, $\Lambda$ by  following the construction in  {cite}`Fournieetal:1999`. 
The first valuation  process tells the marginal impact on future $X$ of a marginal change in one of the initial states analogous to the $\Lambda$ process that we constructed in discrete time. Thus this process has the same number of components as $X$. By initializing the process at one of the alternative coordinate vectors, we again isolate an initial state of interest.[^fn-variational-difference]. 

[^fn-variational-difference]: Our initial condition for $\Lambda_0$ differs from {cite}`Fournieetal:1999` in a superficial way. They treat $\Lambda$ as a matrix with an identity as the initialization. In this way, they consider all of the states of interest simultaneously. We take $\Lambda$ to be a vector and characterize the marginal initial responses one at a time by letting the initial condition be any one of the coordinate vectors.





The drift for the $i^{th}$ component of $\Lambda$ is 

```{math}
\lambda' {\frac {\partial \mu_i} {\partial x} }(x)  
```
and the coefficient on the Brownian increment is 

```{math}
\lambda' \frac {\partial \sigma_i }{\partial x}(x)
```
for $\lambda$ a hypothetical realization of $\Lambda_t$ and $x$ a hypothetical realization of $X_t,$ where $'$ denotes vector or matrix transposition.
The implied evolution of the process $\Lambda^i$ is[^fn-Malliavin]

[^fn-Malliavin]:Since we are working with an instantaneous evolution with Brownian increments, we are implicitly appealing to a formalism known as Malliavin calculus. 

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

Let $\Delta$ be the scalar variational process associated with $Y.$ Then 

```{math}
d \Delta_t = \Lambda_t \cdot \frac {\partial \nu}{\partial x} (X_t)dt + {\Lambda_t}' \frac {\partial \varsigma}{\partial x'} dW_t 
```

Analogous to the discrete-time outcome, the variational dynamics depend explicit on the original diffusion dynamics.  As in discrete time, by initializing the vector $\Lambda_0$ at a coordinate vector, the resulting processes give marginal responses to a corresponding state vector.  

````{prf:example}
:label: ex:continuous
 
Consider the case of linear dynamics:
```{math}
\begin{matrix}
\mu(x)  = {\mathbb A}x  & \sigma(x) = {\mathbb B} \cr
\nu(x)  = {\mathbb D} x &
\varsigma(x)  = {\mathbb F}
\end{matrix} .
```
Then
```{math}
\begin{align*}
\mu^a(x, \lambda) & = \begin{bmatrix} {\mathbb A}x \cr {\mathbb A} \lambda \end{bmatrix} \cr
& \cr
\sigma^a(X, \lambda) & = \begin{bmatrix} {\mathbb B} \cr 0  \end{bmatrix}.
\end{align*} 
 ```
Thus 
```{math}
\Lambda_t = \exp \left( {\mathbb A} t \right) \Lambda_0,
```
and 
```{math} 
\begin{align*}
\Delta_t & = \int_0^t  {\mathbb D} \Lambda_u du + \Delta_0 \cr & = \left[\int_0^t  {\mathbb D} \exp\left( {\mathbb A}u \right) du \right]\Lambda_0 + \Delta_0 \cr & = - {\mathbb A}^{-1} \left[ {\mathbb I} - \exp\left( {\mathbb A} t \right) \right]\Lambda_0 + \Delta_0 .
\end{align*}
```
Given the underlying linearity, the local responses coincide with global responses.  
````

In the calculations that follow, let $\Lambda^j$ be the variational process for which  $\Lambda_0^j$ is a coordinate vector with a one in position $j.$ From the composite processes 
```{math}
\Lambda^a = \begin{bmatrix} \Lambda_1 & ... & \Lambda_n \end{bmatrix}.  
```
Also we will have cause to do a forward shift ${\mathbb S}^\tau $ of these process by which we shift the time units on all of the variables used in the the construction and the initialization period forward $\tau$ time periods.  

 
### Responses to initial shocks

So far, we have characterized stochastic responses to initial changes in the state variables.  From these, we deduce vector of responses to the initial shocks: 

```{math}
\Phi =   \Lambda^a \sigma(X_0)   \hspace{.3cm} \Psi  = \Lambda^a \sigma(X_0)  + \varsigma(X_0) 
```
Under nonlinearity, these responses will be stochastic just as with the state variable perturbations.  
For the special case of linear dynamics given in Example {prf:example}`ex:continuous`, 

```{math}
\begin{align*}
\Phi_t  & =  \exp\left( {\mathbb A} t \right) {\mathbb B}  \cr
\Delta_t  & = - {\mathbb A}^{-1} \left[ {\mathbb I} - \exp\left( {\mathbb A} t \right) \right]{\mathbb B}  + {\mathbb F} ,
\end{align*}
```
which are  the continuous-time counterparts of the familiar impulse responses.  


With Markov diffusions, we also have a state-dependent counterpart to a moving-average representation that is well known from linear time series models.  The resulting formula is known as the Haussmann-Clark-Ocone representation and is given by

```{math}
\begin{align*}
X_t & = \int_0^t   {\mathbb E} \left( {\mathbb S}^u\Phi_{t-u} \mid {\mathfrak F}_u \right) dW_u^j + {\mathbb E} \left( X_t \mid {\mathfrak F}_0 \right) \cr
Y_t & = \int_0^t  {\mathbb E}\left( {\mathbb S}^u \Psi_{t-u}  \mid {\mathfrak F}_u \right) \cdot dW_u^j + {\mathbb E} \left( Y_t \mid {\mathfrak F}_0 \right).
\end{align*} 
```
Note that we form conditional expectations of time shifted  stochastic responses to form the random coefficients in the moving-average representations as given by ${\mathbb E} \left( {\mathbb S}^u\Phi_{t-u} \mid {\mathfrak F}_u \right)$ and ${\mathbb E}\left( {\mathbb S}^u \Psi_{t-u}  \mid {\mathfrak F}_u \right)$.  When the responses turn out not to be stochastic, as in the case of the Remark 2 example, the conditional expectations and the shift are inconsequential.  In this case, we recover the familiar convolution formulas for moving-average representations.   

````{prf:remark}
Many empirical researchers estimate directly what macroeconomists call Jorda projections.  These are implemented by regressing a forward sequence of  a scalar process  on current variable and a shock or particular interest.  One can interpret the ambition as wanting to infer impulse responses from direct regressions of future variables on  the initial ones. One can view the ambition as a way to measure impulse responses.  For instance, the aim could be to infer:
```{math}
{\mathbb E} \left( \Phi_t \mid {\mathcal F}_0 \right)  \hspace{.2cm} {\rm and} \hspace{.2cm}  {\mathbb E} \left( \Psi_t \mid {\mathcal F}_0 \right), \hspace{.2cm} t \ge 0
```
by regressing $X_t$ and $Y_t$ on a measured shock of interest and including additional variables to purge some of the variation in the measured shock.  Many applied papers will include cross terms that are pre-determined in advance of  the shock to accommodate a form of nonlinearity.  For this to be coherent, as our analysis makes clear, one has to think through how the nonlinearity compounds within the stochastic system.  The shock of interest can alter other variables that in turn influence the variable of interest in future time periods.  Our use of variational processes captures this perspective when the ambition is to measure local impacts.
````
