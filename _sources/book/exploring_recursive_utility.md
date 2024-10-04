(chap:recursive)=
# Exploring Recursive Utility 

**Authors:** Borovicka, Jaroslav (NYU), Lars Peter Hansen (University of Chicago) and Thomas J. Sargent (NYU)
$\newcommand{\eqdef}{\stackrel{\text{def}}{=}}$

## Introduction

We explore implications of uncertain-averse economic decision-makers with the recursive utility preference specification of {cite}`KrepsPorteus:1978` and {cite}`EpsteinZin:1989` and counterparts to these preferences that capture concerns about model misspecification as in {cite}`hansensargent:2001` and {cite}`AndersonHansenSargent:2003`.  We investigate the recursrive-utility implications using two different approximation approaches.  First, we consider a continuous-time limiting approximation to a discrete-time preference specification when the underlying shocks are normally distributed.  These results build on a characterization by {cite}`duffieepstein`, but our derivation differs for pedagogical reasons.   We express the representation of the limit in terms of a Brownian motion information structure.  Second, we describe an approximation method that allows macroeconomic uncertainty to have first-order implications. We show how to explore model implications for (nonstandard) first and second-order approximations to equilibria of dynamic stochastic models. We modify first and second order approximations routinely used in the macroeconomics literature in ways that focus on macroeconomic uncertainty.  Our approximations apply to production-based macro-finance models with opportunities to invest in diverse forms of capital. In both approaches, our uncertainty adjustments are captured by a change in probability measure, one that differs from the so-called risk neutral distribution familiar from derivative claims pricing.  We use these frameworks to advance our understanding of alternative preference specifications and their implications for production and asset pricing.





We describe an approximation method that allows macroeconomic uncertainty to have first-order implications. We assume risk averse economic decision makers with recursive utility preferences or closely related preferences that also express concerns about model misspecification. We show how to explore model implications for (nonstandard) first and second-order approximations to equilibria of dynamic stochastic models. We modify first and second order approximations routinely used in the macroeconomics literature in ways that focus on macroeconomic uncertainty. Our uncertainty adjustments are captured by a change in probability measure.  Furthermore, our approximations apply to production-based macro-finance models with opportunities to invest in diverse forms of capital. We use this framework to advance our understanding of alternative preference specifications and their implications for production and asset pricing.




We use an approximation method to explore implications of the recursive utility preference specification of {cite}`KrepsPorteus:1978` and {cite}`EpsteinZin:1989` and counterparts to these preferences that capture concerns about model misspecification. We present formulas for (nonstandard) first and second-order approximations to dynamic, stochastic equilibria for models in which economic agents have such recursive preferences. The approximations build formulations from {cite}`schmitt2004solving` and {cite}`LombardoUhlig:2018`, and we extend them in a way that features the uncertainty contributions more prominently. By design, the implied approximations of stochastic discount factors used to represent market or shadow values reside within the exponential linear quadratic class. This class is known to give tractable formulas for asset valuation over alternative investment horizons. See, for instance, {cite}`AngPiazzesi:2003` and {cite}`borovickahansen14`. Moreover, they are applicable to production-based macro-finance models with investment opportunities in alternative forms of capital.

While these approximations are tractable for the reasons described, for many models researchers may be concerned about some more fundamental aspects of uncertainty that are disguised by these methods.  Indeed, for some models global nonlinearities are more accurately captured by global solution methods.  Nevertheless, the approximations still  provide further understanding of the preferences and their implications for asset pricing and allocation across heterogeneous forms of capital. As a central part of our analysis, we capture the important uncertainty preference contribution as a change in the probability distribution of the underlying economic dynamics. We link this change of measure to the robust preferences specifications of {cite}`hansensargent:2001` and {cite}`AndersonHansenSargent:2003`. The robust preferences formulations build on a robust control literature initiated by {cite}`Jacobson:1973` and {cite}`whittle`.

(sec:approxrecurutil)= 
## Recursive utility valuation process

In this section, we construct  a continuation value  and stochastic discount processes.  These processes are   important constituents of dynamic stochastic models.


### Basic recursion
The homogeneous of degree one representation of recursive utility is 
```{math} 
:label: homog1a
V_t = \left[ (1 - \beta) \left(C_t\right)^{1-\rho} + \beta  \left( R_t \right)^{1-\rho} \right]^{\frac 1 {1-\rho}} 
```
where
```{math} 
:label: homog1b
R_t = \left( {\mathbb E} \left[ \left( V_{t+1} \right)^{1-\gamma} \mid {\mathfrak A}_t \right] \right)^{\frac 1 {1-\gamma}} . 
```
Notice that in equation {eq}`homog1a`, $V_t$ is a homogeneous of degree one function of $C_t$ and $R_t$. In equation {eq}`homog1b`, $R_t$ is a homogeneous of degree one function of another function, namely, $V_{t+1}$ as it varies over date $t+1$ information. 

In equation {eq}`homog1a`, $0 < \beta < 1$ is a subjective discount factor and $\rho$ describes attitudes toward intertemporal substitution. Formally, ${\frac 1 \rho}$ is the elasticity of intertemporal substitution. In equation {eq}`homog1b`, $\gamma$ describes attitudes towards risk.

Continuation values are determined only up to an increasing transformation. For computational and conceptual reasons, we find it advantageous to work with the logarithm ${\widehat V}_t = \log V_t$. The corresponding recursions for ${\widehat  V}_t$ expressed in terms of the logarithm of consumption ${\widehat  C}_t$ are 
```{math} 
:label: value_recur5
{\widehat  V}_t = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t] + \beta \exp \left[(1-\rho) {\widehat  R}_t \right] \right] 
```
where 
```{math} 
:label: value_risk6
{\widehat  R}_t = {\frac 1 {1 - \gamma}} \log {\mathbb E} \left[ \exp \left( (1 - \gamma) {\widehat V}_{t+1} \right) \mid {\mathfrak A}_t \right]. 
```
The right side of recursion {eq}`value_recur5` is the logarithm of a constant elasticity of substitution (CES) function of $\exp({\widehat  C}_t)$ and $\exp({\widehat  R}_t)$.

````{prf:remark}
The limit of ${\widehat R}_t$ as $\gamma$ approaches $1$ is ordinary expected logarithmic utility:
```{math} 
\lim_{\gamma \downarrow 1} \widehat R_t = \lim_{\gamma \downarrow 1} \frac { \log E \left( \exp\left[
(1-\gamma) {\widehat V}_{t+1}\right]
\vert {\mathfrak A}_t \right)}{1-\gamma} = {\mathbb E}\left( \widehat V_{t+1} \vert {\mathfrak A}_t \right) .
```
````

Our approach will be to construct small noise expansions for both ${\widehat V}_t$ and ${\widehat R}_t$ and then to assemble them appropriately. Before doing so, we consider a reinterpretation of {eq}`value_risk6`.

### Preference for robustness

A reinterpretation of the utility recursion and the small-noise expansion approach that we'll deploy comes from recognizing that when $\gamma > 1$, {eq}`value_risk6` emerges from an instance of robust control theory in which $\frac 1 {\gamma - 1}$ is a penalty parameter on entropy relative to alternatives that constrains the alternative probability models that a decision maker considers when evaluating consumption processes. This interpretation originated in work by {cite}`Jacobson:1973` and {cite}`whittle` that was extended and reformulated recursively by {cite}`hansensargentieeetac`.

Let the random variable $N_{t+1} \ge 0$ satisfy ${\mathbb E} \left( N_{t+1} \mid {\mathfrak A}_t \right) = 1$ so that it is a likelihood ratio. Think of replacing the expected continuation value ${\mathbb E} \left( {\widehat V}_{t+1} \mid {\mathfrak A}_t \right)$ by the minimized value of the following problem:
```{math} 
:label: minproblem101
\min_{N_{t+1} \ge 0, {\mathbb E}\left(N_{t+1} \vert {\mathfrak A}_t \right) = 1} {\mathbb E} \left( N_{t+1} {\widehat V}_{t+1} \mid {\mathfrak A}_t \right) + \xi  {\mathbb E} \left( N_{t+1} \log N_{t+1}  \mid {\mathfrak A}_t \right)
```
where $\xi$ is a parameter that penalizes departures of $N_{t+1}$ from unity as measured by relative entropy. Conditional relative entropy for an altered conditional probability induced by applying change of measure $N_{t+1}$ is
```{math} 
{\mathbb E} \left( N_{t+1} \log N_{t+1}  \mid {\mathfrak A}_t \right) \ge 0
```
where, because the function $n \log n$ is convex, the inequality follows from Jensen's inequality. Relative entropy is zero when $N_{t+1} = 1$.  

Relative entropy has a revealing statistical interpretation.  Think of $N_{t+1}$ as the relative conditional likelihood ratio of an alternative model vis-a-vis a baseline model.  Then ${\mathbb E} \left( N_{t+1} \log N_{t+1}  \mid {\mathfrak A}_t \right)$ the expected (conditional) log-likelihood ratio of the alternative model where the expectation is taken using the alternative probability model.  
Small expected log likelihoods are a reflection of the statistical similarity between the two models.  



````{prf:remark}
To solve the minimization problem {eq}`minproblem101`, introduce a Lagrange multiplier, $\ell,$ on the conditional expectation constraint. The Lagrangian problem separates across states, leading to the unconstrained problem:
```{math}
\min_n n v + \xi n \log n + \ell(n - 1)
```
where $n$ is a potential realization of $N_{t+1}$ and $v$ is a realization of $V_{t+1}.$ The first-order conditions are:
```{math}
v + \xi + \xi \log n + \ell = 0.  
```
The solution is
```{math}
n^* = \exp\left[ - \frac 1 \xi \left(v + \ell + \xi \right) \right],
```
and the minimizing objective
```{math}
- \xi \exp\left[ - \frac 1 \xi \left(v + \ell + \xi\right) \right] - \ell .  
```
To complete the solution, we solve for $\ell,$
```{math}
\max_\ell - \xi {\mathbb E} \left( \exp\left[ - \frac 1 \xi \left(V_{t+1} + \ell  + \xi \right) \right] \mid {\mathfrak A}_t \right) - \ell 
```
The first-order conditions are:
```{math}
 {\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right) \exp\left[ - \left(\frac{ \ell + \xi}   \xi\right) \right] -1 = 0.
```
Thus the solution for $\ell$ is
```{math}
\ell^* = \xi \log  {\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right) - 1
```
with a minimized objective given by
```{math}
-  \xi \log  {\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right).
```
The implied minimizer for $N_{t+1}$ is 
```{math} 
:label: Nstar
N_{t+1}^*  = {\frac { \exp \left( - \frac 1 \xi {\widehat  V}_{t+1} \right) }{ {\mathbb E} \left[ \exp \left( - \frac 1 \xi {\widehat  V}_{t+1} \right) \mid {\mathfrak A}_t \right]}} .
```
````

The minimizer of the problem {eq}`minproblem101` given by {eq}`Nstar`
"tilts" probabilities towards low continuation values, what {cite}`bucklew2004` calls a stochastic version of Murphy's law. Notice that the minimized objective satisfies
```{math}
- \xi \log {\mathbb E} \left[ \exp \left( - {\frac 1 \xi} {\widehat  V}_{t+1} \right) \mid
{\mathfrak A}_t \right] = {\widehat R}_t 
```
where ${\widehat  R}_t$ was given previously by equation {eq}`value_risk6` if we set $\xi = {\frac 1 {\gamma - 1}} $.   
It follow from {eq}`Nstar` that 
```{math}
N_{t+1}^* = \exp \left[ - \frac 1 \xi \left({\widehat  V}_{t+1} - {\widehat R}_t \right) \right].
```

The random variable $N_{t+1}^*$ will play a central role in the discussion that follows.  


### Stochastic discount factor process

A stochastic discount factor (SDF) process $S = \{ S_t : t \ge 0 \}$ tells how a consumer responds to small changes in uncertainty and thereby consequently how a consumer values risky payouts. SDF processes have a variety of uses. First, they provide shadow prices that tell how a consumer's uncertainty aversion shapes marginal valuations of risky assets. Second, they shape first-order conditions for optimally choosing financial and physical investments. Third, they underlie tractable formulas for equilibrium asset prices. Fourth, they can help construct Pigouvian taxes for correcting adverse externalities under uncertainty. Fifth, they provide useful tools for assessing effects of small (local) changes in government policies.

To indicate how to deduce an SDF process, we begin by positing that the date zero value of a risky date $t$ consumption payout $\chi_t$ is
```{math}
:label: eqn:price101
\pi_0^t(\chi_t) = E\left[ \left( {\frac {S_t}{S_0}} \right) \chi_t  \Bigr| {\mathfrak A}_0 \right].
```
We compute the ratio $ \frac {S_t}{S_0} $ that appears in formula {eq}`eqn:price101` by evaluating the slope of an indifference curve that runs through both a baseline consumption process $\{C_t\}_{t=0}^\infty$ and a perturbed consumption process
```{math}
\left(C_0 - P_0({\sf q} ), C_1, C_2, \ldots , C_t + {\sf q} \chi_t, C_{t+1}, ... \right).
```
We think of ${\sf q}$ as parameterizing an indifference curve, so $P_0({\sf q})$ expresses how much current period consumption must be reduced to keep a consumer on the same indifference curve after we replace $C_t$ by $C_t + {\sf q} \chi_t$. We set $\pi_0^t(\chi_t)$ defined in equation {eq}`eqn:price101` equal to the slope of that indifference curve:
```{math}
\pi_0^t(\chi_t) = \left. {\frac d {d {\sf q}} } P_0({\sf q}) \right|_{{\sf q} = 0}.
```

The one-period increment in the stochastic discount factor process for recursive utility is:
```{math}
:label: eqn:sdf50
\begin{align*}
\frac {S_{t+1}}{S_t}  = & \beta \left( \frac {C_{t+1}}{C_t} \right)^{-\rho}  \exp\left[ (1-\gamma) \left( {\widehat V}_{t+1} - {\widehat R}_t\right)\right] \exp\left[ (\rho- 1) \left( {\widehat V}_{t+1} - {\widehat R}_t\right)\right]\cr
= & \beta N_{t+1}^* \exp\left({\widehat S}_{t+1} - {\widehat S}_t \right)
\end{align*}
```
where
```{math}
:label: eqn:sdf_change
{\widehat S}_{t+1} - {\widehat S}_t \eqdef  - \rho\left( {\widehat C}_{t+1} + {\widehat C}_t\right)  + (\rho- 1) \left( {\widehat V}_{t+1} - {\widehat R}_{t} \right),
```
and $N_{t+1}^*$ induces the change of probability measure that we described previously as the outcome of robustness problem. (See equation {eq}`Nstar`.) We will use this second formula in {eq}`eqn:sdf50` in what follows.  We think of $\beta$ as a subjective discount factor adjustment, $N_{t+1}^*$ as  a change-of-measure adjustment for uncertainty, and $\exp\left({\widehat S}_{t+1} - {\widehat S}_t \right)$ as an adjustment for the elasticity of intertemporal substitution.  We interpret the transition probality measure induced by $N_{t+1}^*$ as an adjustment for uncertainty in valuation.  
Given the recursive nature of preferences the $t$ period stochastic discount factor increment $\frac {S_t}{S_0}$  as the product of the respective one-period stochastic discount factor increments.  Similarly, we compound the one-period transition uncertainty measure into multiple time-horizon counterparts.




````{prf:remark}
:label: remark_label
To verify formula {eq}`eqn:sdf50`, we compute a one-period intertemporal marginal rate of substitution. Given the valuation recursions {eq}`value_recur5` and {eq}`value_risk6`, we construct two marginal utilities familiar from CES and exponential utility:
```{math}
mc = (1 - \beta)\left(   c   \right)^{-\rho} \exp\left[ (\rho - 1) {\hat v} \right]
```
```{math}
m{\hat r} = \beta  \exp[(1-\rho) \left({\hat r} - {\hat v}   \right) ]
```
From the certainty equivalent formula, we construct the marginal utility of the next-period logarithm of the continuation value:
```{math}
m{\hat v}^+=  \exp \left[(1-\gamma ) \left({\hat v}^+ -  {\hat r} \right) \right]
```
where the $+$ superscript is used to denote the next-period counterpart. In addition, the next-period marginal utility of consumption is
```{math}
mc^+ = (1 - \beta)\left(   c ^+  \right)^{-\rho} \exp\left[ (\rho - 1) {\hat v}^+\right] 
```
Putting these four formulas together using the chain rule for differentiation gives a marginal rate of substitution:
```{math}
\frac {(mr) (mv^+) (mc^+)}{mc} = \beta \left( \frac {c^+}{c} \right)^{-\rho} \exp \left[(1-\gamma ) \left({\hat v}^+ -  {\hat r} \right) \right]
\exp \left[(\rho - 1) \left({\hat v}^+ -  {\hat r} \right) \right]. 
```
Now let ${\hat v}^+ =  {\widehat V}_{t+1}$, $c^+ = C_{t+1}$, $C_t = c$ and ${\hat r} = {\widehat R}_t$ to obtain the formula for the one-period stochastic discount factor {eq}`eqn:sdf50`.
````

## Continuous-time limit

We now  explore a continuous-time limit that approximates a discrete-time specification. Since we continue to work with normal shocks, in the continuous-time counterpart to these shocks are Brownian increments.
The continuation value in continuous time will evolve as:
```{math}
d {V}_t = V_t {\mu}_t^V dt  + V_t\sigma_t^V \cdot dW_t
```
for some drift (local mean), $V_t \mu_t^V$ and some local shock exposure vector $V_t \sigma_t^V$, where $\{W_t  : t \ge 0\}$ is a multivariate Brownian motion.  The scaling for the local evolution coefficients by $V_t$ is done for convenience, where the continuation value process is presumed to be positive.  As in discrete time, it is convenient to work with the logarithm of the continuation value process (a strictly increasing transformation.  The implied evolution is  
```{math}
d {\widehat V}_t = {\hat \mu}_t^V dt  + \sigma_t^V \cdot dW_t
 ``` 
where ${\hat \mu}_t = {\mu}_t^V - {\frac 1 2} \mid \sigma_t^V \mid^2.$ This adjustment follows from the well known Ito formula.  

### Discrete-time approximation 

To study the utility recursion, start with a discrete-time specification:

```{math}
\begin{align*}
 & \frac 1 {1-\rho} \log 
 \left[(1 - \beta_\epsilon) \exp\left[(1-\rho)\left({\widehat C}_t - {\widehat V}_t \right)\right] + \beta_\epsilon \exp \left[(1-\rho)\left({\widehat R}_t -{\widehat V}_t \right) \right] \right] = 0 \cr
& {\widehat  R}_t - {\widehat V}_t =  \frac 1 {1 - \gamma}  \log {\mathbb E}  \left( \exp\left[ (1-\gamma)\left( {\widehat V}_{t+\epsilon} - {\widehat V}_t  \right)\right] \mid {\mathfrak A}_{t}\right) 
\end{align*}
```
where $\beta_\epsilon = \exp(- \delta \epsilon)$ and $\delta>0$ is instantaneous subjective rate of discount. Consider the time derivative of the second recursion:
```{math}
\begin{align*} 
\frac d {d\epsilon} \left. \left({\widehat  R}_t - {\widehat V}_t\right) \right|_{\epsilon = 0} & =  \frac d {d\epsilon} \left. \frac 1 {1 - \gamma}  \log {\mathbb E}  \left( \exp\left[ (1-\gamma)\left( {\widehat V}_{t+\epsilon} - {\widehat V}_t  \right)\right] \mid {\mathfrak A}_{t}\right)\right|_{\epsilon = 0} \cr
& = {\hat \mu_t^V} 
  + \frac {(1 - \gamma)} 2  |\sigma_t^V|^2 \cr
 & =  \mu_t^V 
  - \frac \gamma 2 |\sigma_t^V|^2 .
\end{align*}
```
by local log normality.  We turn now to the first recursion and compute time derivatives in three steps.  First, we 
evaluate the term inside the logarithm as $\epsilon$ tends to zero :
```{math}
\lim_{\epsilon \downarrow 0} \exp\left[ (1-\rho)\left({\widehat  R}_t - {\widehat V}_t\right) \right] = 1. 
```
This term is in the denominator  as implied by the derivative   with respect to a logarithm.  Second, we 
differentiate the term inside the logarithm with respect to $\epsilon$ as contributed by $\beta_\epsilon$:
```{math}
\delta \exp\left[(1-\rho)\left({\widehat C}_t - {\widehat V}_t \right)\right] - \delta .
```
Third, we differentiate the term inside the logarithm with respect to $\epsilon$ as contributed by 
```{math}
\frac d {d\epsilon} \exp \left[(1-\rho)\left({\widehat R}_t -{\widehat V}_t \right) \right] = {\hat \mu}_t^V + \frac {(1-\gamma)} 2 |\sigma_t^V|^2.
```
Putting together the derivative components together gives:
```{math}
:label: cont_recursion
\begin{align*}
 0 = \hspace{.2cm} &\frac {\delta 
\left[\left(\frac {C_t}{V_t} \right)^{1-\rho} - 1\right]}{1 - \rho} + {\hat \mu}_t^V + \frac {(1-\gamma)} 2 |\sigma_t^V|^2\cr
= \hspace{.2cm} &\frac {\delta 
\left[\left(\frac {C_t}{V_t} \right)^{1-\rho} - 1\right]}{1 - \rho} + {\mu}_t^V - \frac {\gamma} 2 |\sigma_t^V|^2.
\end{align*}
```

Notice that this relation imposes a restriction across the local mean, ${\mu}_t^V$ and local variance, $ |\sigma_t^V|^2$ of the  continuation value. {cite}`duffieepstein` refer to $\gamma$ as a  *variance multiplier* where larger values of $\gamma$ imply a more substantial adjustment for local volatility.  As in discrete-time, the $\rho = 1$ is interesting special case represented as;
```{math}
 0 = \delta \left( {\widehat V}_t - {\widehat C}_t \right)  + {\hat \mu}_t^V + \frac {(1-\gamma)} 2 |\sigma_t^V|^2.
```

### Robustness to misspecification

To investigate an aversion to model misspecification in continuous time, 
we now treat the distribution of $\{W_t : t\ge0\}$ as uncertain.  We allow for  probability measures that entertain possible Brownian motions with local means or drifts that are history dependent.  



We start by considering positive martingales $\{ M^H_t : t \ge 0\}$ parameterized by 
alternative $\{ H_t : t \ge 0\}$ processes with the same dimension as the underlying Brownian motion.  
The martingales have local evolutions:
```{math}
d M_t^H = M_t  H_t \cdot dW_t,
```
and we initialize them at $M_0^H = 1.$.  Observe  that by applying Ito's formula,   $\log M$ evolves as:
```{math}
d \log M_t^H = -{\frac 1 2} \left| H_t \right|^2 dt + H_t \cdot dW_t.
```
We use these martingales as relative densities or likelihood ratios.  


Write the discrete-time counterpart as
```{math}
\log M_{t+\epsilon}^H - \log M_t^H  =  -{\frac \epsilon 2} \left| H_t \right|^2  + H_t \cdot \left(W_{t+\epsilon} -W_t\right) 
```
Let $w$ be a realized $W_{t+\epsilon} - W_t$ and $h$ be a realization $H_t$.  Then $\log M_{t+\epsilon} - \log M_t$ contributes
$  -{\frac \epsilon 2}h'h + h\cdot w $ to the log-likelihood.  The standard normal density for $\sqrt{\epsilon } \left(W_{t+\epsilon} -W_t\right)$ contributes $-{\frac 1 \epsilon} w'w - \log (2 \pi \epsilon)$.  Put there two components together, we have a log-likelihood:
```{math}
-{\frac \epsilon 2}h'h + h\cdot w - \frac 1 {2 \epsilon} w'w - \log (2 \pi \epsilon) = - \frac 1 {2 \epsilon} (w - \epsilon h)' (w - \epsilon h) -  \log (2 \pi \epsilon).
```
The altered conditional density has mean $\epsilon h$, which is the realized value of $\epsilon H_t$ with the same conditional covariance matrix as before. Moreover, the conditional expectation of 
```{math}
{\mathbb E} \left[\left(\frac {M_{t+\epsilon}^H}{M_t^H}\right)\left(\log M_{t+\epsilon}^H - \log M_t^H \right) \mid {\mathfrak A}_t \right] = {\frac \epsilon 2} \left| H_t \right|^2, 
```
which measures the statistical divergence or relative entropy between original and altered conditional probabilities.


For the continuous-time limit,  
under the $H$ change of probability measure:
```{math}
dW_t = H_t dt + d{\widetilde W}_t^H
```
where ${\widetilde W}^H$ is a standard Brownian motion.  Thus the potential changes of probability measures induce a \alert{local means} or drift $H$ processes to  the Brownian motion.  The continuous-time counterpart to conditional relative entropy at time $t$ is ${\frac 1 2} \left| H_t \right|^2.$

We may justify focusing on drift distortions for Brownian increments because or our imposition of absolute continuity of the alternative probabilities with with respect to the baseline specification of a multivariate standard Brownian motion.  This is an implication of the  Girsanov Theorem.  

We are now in a position to deduce a robustness adjustment in continuous time.
Consider formula {eq}`cont_recursion` when $\gamma = 1$ modified for a potential change in the probability measure
```{math}
0 = \frac {\delta 
\left[\left(\frac {C_t}{V_t} \right)^{1-\rho} - 1\right]}{1 - \rho} + {\hat \mu}_t^V + \sigma_T^V \cdot H_t 
```
Modify this equation to include minimization over $H_t$ subject to a relative entropy penalty $\frac \xi 2 |H_t|^2:$
```{math} 
0  = \min_{H_t} \frac {\delta 
\left[\left(\frac {C_t}{V_t} \right)^{1-\rho} - 1\right]}{1 - \rho} + {\hat \mu}_t^V + \sigma_t^V\cdot  H_t + \frac {{\xi}} 2 H_t \cdot H_t.  
```
The minimizing solution is 
```{math}
H_t^* =  - \frac 1 {\xi} \sigma_t^V 
```
with a minimized objective 
```{math}
0 = \frac {\delta 
\left[\left(\frac {C_t}{V_t} \right)^{1-\rho} - 1\right]}{1 - \rho} + {\hat \mu}_t^V   - \frac 1 {2 \xi} \left| \sigma_t^V \right|^2
```
Notice that this  agrees with formula {eq}`cont_recursion` for $\gamma - 1= 1/\xi$.  The explicit link is entirely consistent with our discrete-time equivalence result.  By taking the continuous-time limit, we are able to focus our misspecification analysis on changing local means of the underlying Brownian increments.  
 
### Uncertainty pricing 

For the purposes of valuation, we compound the equilibrium version of $\{ H_t^* : t \ge 0\}$ to give a  exponential  martingale:
```{math}
M_t^* = \exp \left( \int_0^t H_\tau^* dW_{d\tau}  - {\frac 1 2} \int_0^t \left| H_\tau^* \right|^2 d\tau \right)
```
provided that the constructed process is a martingale.[^local]  With this construction we interpret $-H_t^*$ as the vector of local uncertainty prices that give compensations for exposure to Brownian increment uncertainty.  These compensations are expressed as changes in conditional means under the baseline distribution as is typical in continuous-time asset pricing. 

[^local]: In general, this exponential martingale formula produces a local martingale with conditional expectations that might decline over time.  There are a variety of sufficient conditions that may be checked to verify that the constructed process is actually a martingale with unit expectation.



##  Small noise expansion of dynamic stochastic equilibria

We next consider a different type of characterization  that sometimes gives good approximations for dynamic stochastic equilibrium models.  While the approximations build  from derivations in {cite}`schmitt2004solving` and {cite}`LombardoUhlig:2018`, and we extend them in a way that features the uncertainty contributions more prominently and are reflected even in first-order contributions.  By design, the implied approximations of stochastic discount factors used to represent market or shadow values reside within the exponential linear quadratic class. This class is known to give tractable formulas for asset valuation over alternative investment horizons. See, for instance, {cite}`AngPiazzesi:2003` and {cite}`borovickahansen14`. Moreover, they are applicable to production-based macro-finance models with investment opportunities in alternative forms of capital.

While these approximations are tractable for the reasons described,  researchers may be concerned about some more fundamental aspects of uncertainty that are disguised by these methods.  Indeed, for some models  nonlinearities are more accurately captured by global solution methods.  Nevertheless, the approximations still  provide further understanding of the preferences and their implications for asset pricing in endowment and production economies. 

### Approximate state dynamics
We follow {cite}`LombardoUhlig:2018` by considering the following class of stochastic processes indexed by a scalar perturbation parameter $\mathsf{q}$:[^footnote_LombardoUhlig]

```{math}
:label: eq:xlom
X_{t+1}\left( \mathsf{q}\right) =\psi \left[ X_{t}\left( \mathsf{q}\right) ,\mathsf{q}W_{t+1},\mathsf{q}\right] .
```
Here $X$ is an $n$-dimensional stochastic process and $\{W_{t+1}\}$ is an i.i.d.~normally distributed random vector with conditional mean vector $0$ and conditional covariance matrix $I$. We parameterize this family so that ${\sf q} = 1$ gives the model of interest.

We denote a zero-order expansion ${\sf q} = 0$ limit as:
```{math}
:label: eqn:exporder0
X^{0}_{t+1}  = \psi \left( X_{t}^0 ,0,0\right),
```
and assume that there exists a second-order expansion of $X_{t}$ around $\mathsf{q} = 0$:
```{math}
:label: eqn:exporder2
X_{t} \approx X_{t}^0+\mathsf{q}X_{t}^1+\frac{{\sf q} ^2}{2}X_{t}^2 
```
where $X_t^1$ is a first-order contribution and $X_t^2$ is a second-order contribution. In other words, the stochastic processes $X^j$, $j=0,1,2$ are appropriate derivatives of $X$ with respect to the perturbation parameter ${\sf q}$ evaluated at ${\sf q} = 0.$

In the remainder of this chapter, we shall construct instances of the second-order expansion {eq}`eqn:exporder2` in which the generic random variable $X_t$ is replaced, for example, by the logarithm of consumption, a value function, and so on. I

Processes $X_t^j, j=0, 1, 2$ have a recursive structure: first compute the stochastic process $X_{t}^0,$  then the process $X_{t}^1$ next (it depends on $X_t^0$), and finally the process $X_{t}^2$ (it depends on both $X_{t}^0$ and $X_{t}^1$).

In this chapter, we use a prime $'$ to denote a transpose of a matrix or vector. When we include $x'$ in a partial derivative of a scalar function it means that the partial derivative is a row vector. Consistent with this convention, let $\psi^i_{x'}$, the $i^{th}$ entry of $\psi_{x'}$, denote the row vector of first derivatives with respect to the vector $x$, and similarly for $\psi^i_{w'}$. Since ${\sf q}$ is scalar, $\psi_{\sf q}^i$ is the scalar derivative with respect to ${\sf q}$. Derivatives are evaluated at $X_t^0$, which in many examples is invariant over time, unless otherwise stated. This invariance follows when we impose a steady state on the deterministic system.

The first-derivative process obeys a recursion
```{math}
:label: exprecur1
X_{t+1}^1 =  \begin{bmatrix} \psi_{x'}^1 \\ \psi_{x'}^2 \\ \vdots \\ \psi_{x'} ^n \end{bmatrix} X_t^1 + \begin{bmatrix} \psi_{w'}^1 \\ \psi^2_{w'} \\ \vdots \\ \psi^n_{w'} \end{bmatrix} W_{t+1} + \begin{bmatrix} \psi^1_{\sf q}  \\ \psi^2_{\sf q}  \\ \vdots \\ \psi^n_{\sf q}  \end{bmatrix}
```
that we can write compactly as the following *a first-order vector autoregression*:
```{math}
X_{t+1}^1 = \psi _{x'}X_{t}^1+\psi _{w'}W_{t+1}+\psi_{\sf q}
```
We assume that the matrix $\psi_x'$ is stable in the sense that all of its eigenvalues are strictly less than one in modulus.

It is natural for us to denote second derivative processes with double subscripts. For instance, for the double script used in conjunction with the second derivative matrix of $\psi^i$, the first subscript without a prime ($'$) reports the row location; second subscript with a prime ($'$) reports the column location. Differentiating recursion {eq}`exprecur1` gives:

```{math}
:label: exprecur2
\begin{split}
X_{t+1}^2 =  \hspace{.2cm} \psi_{x'} X_t^2 &+ \begin{bmatrix} X_t^{1'} \psi_{xx'}^1 X_t^1 \\ X_t^{1'} \psi_{xx'}^2  X_t^1\\ \vdots \\ X_t^{1'} \psi_{xx'}^n X_t^1 \end{bmatrix} +  2 \begin{bmatrix} X_t^{1'} \psi_{xw'}^1 W_{t+1}  \\ X_t^{1'} \psi_{xw'}^2  W_{t+1} \\ \vdots \\ X_t^{1'} \psi_{xw'}^n W_{t+1}  \end{bmatrix} + \begin{bmatrix} {W_{t+1}}' \psi_{ww'}^1 W_{t+1}  \\ {W_{t+1}}' \psi_{ww'}^2  W_{t+1} \\ \vdots \\ {W_{t+1}}' \psi_{ww'}^n W_{t+1}  \end{bmatrix} 
 \\ &+  2 \begin{bmatrix} \psi^1_{{\sf q}  x'}X_{t}^1 \\  \psi^2_{{\sf q}  x'}X_{t}^1 \\ \vdots \\ \psi^n_{{\sf q}  x'} X_t^1 \end{bmatrix}
+ 2 \begin{bmatrix} \psi_{{\sf q}  w' }^1W_{t+1} \\ \psi_{{\sf q}  w' }^2 W_{t+1} \\ \vdots \\ \psi^n_{{\sf q}  w'} W_{t+1} \end{bmatrix}
+ \begin{bmatrix} \psi^1_{{\sf q}  {\sf q} } \\ \psi^2_{{\sf q} {\sf q} } \\ \vdots \\ \psi^n_{{\sf q} {\sf q} } \end{bmatrix}.
\end{split}
```
Recursions {eq}`exprecur1` and {eq}`exprecur2` have a linear structure with some notable properties. The law of motion for $X^0$ is deterministic and is time invariant if {eq}`eq:xlom` comes from a stationary $\{X_t\}$ process. The dynamics for $X^2$ are nonlinear only in $X^1$ and $W_{t+1}$. Thus, the stable dynamics for $X^1$ that prevail when $\psi_x$ is a stable matrix imply stable dynamics for $X^2$. 

[^footnote_LombardoUhlig]: {cite}`LombardoUhlig:2018` provides a discussion of how their approach builds on more general perturbation methods as discussed by {cite}`Holmes:2012` and {cite}`Judd:1998`.

````{prf:remark}
Perturbation methods have been applied to many rational expectations models in which partial derivatives of $ \psi $ with respect to $\mathsf{q}$ are often zero.[^footnote1] However, derivatives of $ \psi $ with respect to $\mathsf{q}$ aren't zero in production-based equlibrium models with the robust or recursive utility specifications that we shall study here.
````

Let $C$ denote consumption and ${\widehat C}$ the logarithm of consumption. Suppose that the logarithm of consumption evolves as:

```{math}
{\widehat C}_{t+1} - {\widehat C}_t = \kappa (X_t, {\sf q} W_{t+1}, {\sf q}).
```

Approximate this process by:

```{math}
:label: eqn:exporder2c
{\widehat C}_{t+1} - {\widehat C}_t \approx 
{\widehat C}_{t+1}^0 - {\widehat C}_t^0 + {\sf{q}} \left({\widehat C}_{t+1}^1 - {\widehat C}_t^1 \right)
+ \frac{{\sf q}^2} 2 \left({\widehat C}_{t+1}^2 - {\widehat C}_t^2\right) 
```

where 

```{math}
{\widehat C}_{t+1}^0  - {\widehat C}_t^0  = \hspace{.2cm}  \kappa (X_t^0, 0, 0) := \eta_0^c 
```

```{math}
{\widehat C}_{t+1}^1  - {\widehat C}_t^1  = \hspace{.2cm} \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1} + \kappa_q 
```

```{math}
\begin{align*}
{\widehat C}_{t+1}^2  - {\widehat C}_t^2  =  &  \hspace{.2cm}  \kappa_{x'} X_t^2 + {X_{t}^1}' \kappa_{x,x'} X_t^1 
+ 2 \kappa_{q,x'} X_t^1 + \kappa_{qq} \cr 
& + 2
{X_{t}^1}' \kappa_{xw'} W_{t+1} + {W_{t+1}}' \kappa_{ww'} W_{t+1} 
+  +  2 \kappa_{qw'} W_{t+1}  .
\end{align*}
```

In models with endogenous investment and savings, the consumption dynamics and some of the state dynamics will emerge as the solution to a dynamic stochastic equilibrium model. We use the approximating processes {eq}`eqn:exporder2` and {eq}`eqn:exporder2c` as inputs into the construction of an approximating continuation value process and its risk-adjusted counterpart for recursive utility preferences.

[^footnote1]: See, for instance, {cite}`schmitt2004solving`.


## Incorporating preferences with enhanced uncertainty concerns

To approximate the recursive utility process, we deviate from common practice in macroeconomics by letting the risk aversion or robust parameter in preferences depend on ${\sf q}:$
```{math}
\xi = {\sf q} \xi_o    \hspace{1 cm}  \gamma - 1 =  \frac {\gamma_o - 1} {\sf q} 
```
The aversion to model misspecification or the aversion to risk moves inversely with the parameter ${\sf q}$ when we embed the model of interest within a parameterized family of models. In effect, the variable ${\sf q}$ is doing double duty. Reducing ${\sf q} > 0$ limits the overall exposure of the economy to the underlying shocks. This is offset by letting the preferences include a greater aversion to uncertainty. This choice of any expansion protocol has significant and enlightening consequences for continuation value processes and for the minimizing $N$ process used to alter expectations. It has antecedents in the control theory literature, and it has the virtue that implied uncertainty adjustments occur more prominently at lower-order terms in the approximation.  

### Order-zero

Write the order-zero expansion of {eq}`value_recur5` as
```{math}
{\widehat  V}_t^0 = {\frac 1 {1 - \rho}} \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t^0] + \beta \exp \left[(1-\rho) {\widehat  R}_t^0 \right] \right]
```
```{math}
{\widehat  R}_t^0  =  {\widehat  V}_{t+1}^0,
```
where the second equation follows from noting that randomness vanishes in the limit as $\mathsf{q}$ approaches $0$.

For order zero, write the consumption growth-rate process as
```{math}
{\widehat C}_{t+1}^0 - {\widehat  C}_t^0  = \eta_c^0 .
```
The order-zero approximation of {eq}`value_recur5` is:
```{math}
{\widehat  V}_t^0 - {\widehat  C}_t^0 = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta)  + \beta \exp \left[(1-\rho) \left( {\widehat  V}_{t+1}^0 - {\widehat C}_{t+1}^0 + \eta_c^0 \right)  \right] \right]
```
We guess that ${\widehat  V}_t^0 - {\widehat  C}_t^0 = \eta_{v-c}^0$ and will have verified the guess once we solve:
```{math}
\exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right] = (1 - \beta) + \beta \exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right]\exp \left[ (1 - \rho) \eta_c^0 \right].
```
This equation implies
```{math}
:label: zero_formula
\exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right] = {\frac {1 - \beta} { 1 - \beta \exp \left[ (1 - \rho) \eta_c^0 \right]}} .
```
Equation {eq}`zero_formula`
 determines $ \eta_{v - c}^0$ as a function of $\eta_c^0$ and the preference parameters $\rho, \beta$, but not the risk aversion parameter  $\gamma$ or its robust counterpart $\xi$. Specifically,
```{math}
:label: zero_formula2
 \eta_{v-c}^0 = \frac {\log (1 - \beta)  - \log \left(1 - \beta \exp \left[ (1 - \rho) \eta_c^0 \right]\right)}{1 - \rho} 
 ```
In the limiting $\rho =1$ case, 
```{math}
\eta_{v-c}^0 = \left(\frac \beta {1 - \beta} \right) \eta_c^0.
```

### Order-one

We temporarily take ${\widehat  R}_t^1 - {\widehat C}_t^1$ as given.  We construct a first-order approximation to the nonlinear utility recursion {eq}`value_recur5`
```{math}
:label: eqn:lambda1
{\widehat  V}_t^1 - {\widehat C}_t^1 = \lambda \left({\widehat  R}_t^1 - {\widehat C}_t^1\right)
```
where
```{math}
\begin{align*}
\lambda  = & \left[ {\frac {\beta \exp \left[(1-\rho) \left(\eta_{ v - c}^0 +\eta_c^0 \right ) \right]} { (1 - \beta)  + \beta \exp \left[(1-\rho) \left(\eta_{ v - c}^0 + \eta_c^0  \right ) \right]}} \right] \cr \cr
= & \left[ {\frac {\beta \exp \left[(1-\rho) \eta_c^0\right]} { (1 - \beta) \exp\left[- (1-\rho) \eta_{ v - c}^0\right]   + \beta \exp \left[(1-\rho) \eta_c^0 \right]}} \right] \cr \cr
= & \left[ {\frac {\beta \exp \left[(1-\rho) \eta_c^0 \right]} {  1 - \beta \exp \left[ (1 - \rho)\eta_c^0 \right]  + \beta \exp \left[(1-\rho) \eta_c^0 \right]}} \right] \cr \cr
= & \beta \exp \left[(1-\rho) \eta_c^0 \right].
\end{align*}
```
Notice how the parameter $\rho$ influences the weight $\lambda$ when $\eta_c \neq 0$, in which case the log consumption process displays growth or decay. The condition $\lambda  <1$ restricts the subjective discount rate, $-\log \beta,$ relative to the consumption growth rate $\eta_c$ since
```{math}
(1- \rho) \eta_c^0 <  - \log \beta .
```
When $\rho < 1$, the subjective discount rate has a positive lower bound in contrast to the case in which $\rho \ge 1$.  


To facilitate computing some useful limits we construct:
```{math}
:label: tilde_transform
{\widetilde V}_t \eqdef {\frac {{\widehat  V}_t - {{\widehat  V}_t^0}}{{\sf q} }} 
```
```{math}
:label: tilde_R
{\widetilde R}_t \eqdef  {\frac {{\widehat  R}_t - {{\widehat  V}_{t+1}^0}}{{\sf q} }},
```
which we assume remain well defined as ${\sf q} $ declines to zero, with limits denoted by ${\widetilde V}_t^0$ and ${\widetilde  R}_t^0$.
Importantly, 
```{math}
{\widetilde R}_t =  \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right] \mid {\mathfrak A}_t \right),
```
Taking limits as ${\sf q}$ declines to zero:
```{math} 
{\widehat R}_t^1 =  \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widehat V}_{t+1}^1 \right] \mid {\mathfrak A}_t \right)
```
Subtracting ${\widehat C}_t^1$ from both sides gives:
```{math}
:label: first_order_risk
\begin{align}
& {\widehat R}_t^1 - {\widehat C}_t^1  =  \cr  & \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o) \left[ \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1} ^1 \right) +   \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t} ^1 \right)\right]  \right] \mid {\mathfrak A}_t \right).
\end{align}
```
Substituting formula {eq}`first_order_risk` into the right side of {eq}`eqn:lambda1` gives the recursion for the first-order continuation value:
```{math}
:label: first_recursive_update
\begin{align}
& {\widehat V}_t^1 - {\widehat C}_t^1  =   \cr 
& \left( {\frac \lambda {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o) \left[ \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1} ^1 \right) +   \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t}^1 \right) \right] \right] \mid {\mathfrak A}_t \right).
\end{align}
```

````{prf:remark}
:label: remark_solve_first  
We produce a solution by "guess and verify." Suppose that 
```{math} 
:label: first-solution 
{\widehat V}_t^1 - {\widehat C}_t^1 = {\upsilon_1}' X_t^1 + \upsilon_0. 
```
It follows from {eq}`first_recursive_update` that
```{math} 
:label: first_formula 
{\upsilon_1} = \lambda \left(\left(\psi_{x'}\right)'{\upsilon_1} + \kappa_{x} \right), 
\upsilon_0 = \lambda \left( \upsilon_0 + {\upsilon_1}' \psi_{q} + \kappa_q + {\frac {(1 - \gamma_0)} 2} 
\left| {\upsilon_1}' \psi_{w'} + \kappa_{w'} \right|^2 \right). 
```
Deduce the second equation by observing that 
```{math}
 \exp \left[ (1 - \gamma_0)  \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1}^1 \right) + (1-\gamma_0)  \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t}^1 \right) \right]
 ```
 is distributed as a log normal. The solutions to equations {eq}`first_formula` are:
```{math}
\upsilon_1 = \lambda  \left( I - \lambda \psi_{x'} \right)^{-1} \kappa_{x'}, 
\upsilon_0 = {\frac \lambda {(1 - \lambda)}} \left({\upsilon_1}'  \psi_{q} +  \kappa_{q} \right) +  
\frac {\lambda(1 - \gamma_0)} {2(1 - \lambda)} 
\left| {\upsilon_1}' \psi_{w'} + \kappa_{w'}  \right|^2. 
```
The continuation value has two components. The first is:
```{math} 
{\upsilon_1}'X_t^1  + {\frac \lambda {(1 - \lambda)}} \left({\upsilon_1}'  \psi_{q} +  \kappa_{q} \right)
= {\mathbb E} \left[ \sum_{j=1}^\infty \lambda^j \left({\widehat C}_{t+j}^1 - {\widehat C}_{t+j-1} ^1\right) \mid {\mathfrak A}_t  \right] 
```
and the second component is a constant long-run risk adjustment given by:
```{math} 
\frac {\lambda(1 - \gamma_0)} {2(1 - \lambda)} 
\left| {\upsilon_1}' \psi_{w'} + \kappa_{w'} \right|^2.
```
This second term is the variance of 
```{math} 
:label: discount_future 
 {\mathbb E} \left[\sum_{j=1}^\infty \lambda^j \left({\widehat C}_{t+j }^1 - {\widehat C}_{t+j-1} ^1\right) \mid {\mathfrak A}_{t+1}  \right] = (1-\lambda) {\mathbb E} \left[\sum_{j=1}^\infty \lambda^j \left({\widehat C}_{t+j }^1 - {\widehat C}_{t} ^1\right) \mid {\mathfrak A}_{t+1}  \right] 
```
conditioned on ${\mathfrak A}_t$ scaled by $\frac {\lambda(1 - \gamma_0)} {2(1 - \lambda)}$. 
````

````{prf:remark}
:label: remark:zero
The formula for $ \upsilon_1 $ depends on the parameter $ \rho $. Moreover, $ \upsilon_1 $ has a well-defined limit as $ \lambda $ tends to unity as does the variance of {eq}`discount_future`. This limiting variance:
```{math}
\lim_{\lambda \rightarrow 1} \left| {\upsilon_1}' \psi_{w'} + \kappa_{w'} \right|^2.
```
converges to the variance of the martingale increment of ${\widehat C}^1$.  
````

````{prf:remark}
:label: remark:first
Consider the logarithm of the uncertainty-adjusted continuation value approximated to the first order.  Note that from {eq}`first-solution`,
```{math}
{\widehat V}_{t+1}^1 - {\widehat C}_t^1  = {\upsilon_1}' X_{t+1}^1 + \upsilon_0 + \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1}.
```
Substitute this expression into formula {eq}`first_order_risk` and use the formula for the mean of random variable distributed as a log normal to show that
```{math}
{\widehat V}_{t+1}^1 - {\widehat R}_t^1 = \left( {\upsilon_1}'\psi{_w'} + \kappa_{w'} \right) W_{t+1} - \left( \frac{1 - \gamma_o}{2} \right)
\vert  {\upsilon_1}'\psi{_w'} + \kappa_{w'} \vert^2.
```
````



Associated with the first-order approximation, we construct:
```{math}
N_{t+1}^0  \eqdef {\frac{ \exp \left[ (1 - \gamma_o)   {\widetilde V}_{t+1}^0 \right]}{{\mathbb E}\left( \exp \left[ (1 - \gamma_o)   {\widetilde V}_{t+1}^0 \right] \mid {\mathfrak A}_t \right)}} = \frac{ \exp \left[ (1 - \gamma_o)   {\widehat V}_{t+1}^1 \right]}{{\mathbb E} \left( \exp \left[ (1 - \gamma_o)   {\widehat V}_{t+1}^1 \right] \mid {\mathfrak A}_t \right)}.
```



Equation {eq}`first_order_risk` is a standard risk-sensitive recursion applied to log-linear dynamics. For instance, see {cite}`tallarini`'s paper on risk-sensitive business cycles and {cite}`HansenHeatonLi:2008`'s paper on measurement and inference challenges created by the presence of long-term risk.[^footnoteHHL] Both of those papers assumed a logarithmic one-period utility function, so that for them $\rho=1.$ Here we have instead obtained the recursion as a first-order approximation without necessarily assuming log utility. Allowing for $\rho$ to be different than one shows up in both the order zero and order one approximations, as reflected in {eq}`zero_formula2` and {eq}`first_recursive_update`, respectively. In accordance to  {eq}`first_recursive_update`, for the first-order approximation the parameter $\lambda = \beta$ when $\rho = 1$. But otherwise, it is different. Equation {eq}`first_order_risk` also is very similar to a first-order approximation proposed in {cite}`RestoyWeil:2011`. Like formula {eq}`first_order_risk`, {cite}`RestoyWeil:2011` allow for $\rho \ne 1.$ In contrast, our equation has an explicit constant term coming from the uncertainty adjustment, and we have explicit formula for $\lambda$ that depends on preference parameters and the consumption growth rate.



[^footnoteHHL]:Consistent with overall message of the paper, the {cite}`HansenHeatonLi:2008` predictability evidence turned out to be "fragile" and was modified and updated in {cite}`HansenSargent:2021` Appendix B. This same appendix suggests a way to deduce a statistical approximation to the first order dynamics of {cite}`bansalyaron2004` from a more general VAR representation of the consumption dynamics.








````{prf:remark}
:label: remark:second
The calculation reported in {prf:ref}`remark:first` implies that
```{math}
\begin{align}
\log N_{t+1}^0 = (1-\gamma_o)\left({\widehat V}_{t+1}^1- {\widehat R}_t^1\right) =  \hspace{.2cm} & (1- \gamma_o) \left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right) W_{t+1} \cr \hspace{.2cm} & - \frac{(1 - \gamma_o)^2}{2} \vert {\upsilon_1}'\psi_{w'} + \kappa_{w'} \vert^2
\end{align}
```
As a consequence, under the change in probability measure induced by $N_{t+1}^0,$ $W_{t+1}$ has a mean given by
```{math}
\begin{align}
\mu^0 \eqdef \hspace{.2cm} (1 - \gamma_o)  \left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right)' \cr = \hspace{.2cm}\left( \frac 1 {\xi_o} \right) \left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right)'
\end{align}
```
and with the same covariance matrix given by the identity. This is an approximation to robustness adjustment expressed as an altered distribution of the underlying shocks. It depends on $\gamma_o - 1 = {\frac{1}{\xi_o}}$ as well as the state dynamics as reflected by $\upsilon_1$ and by the shock exposure vectors $\psi_{w'}$ and $\kappa_{w'}$.   As we will see, this change of measure plays a role in the higher-order approximation, but it also gives a low-order representation of the implied shadow or market one-period compensation for exposure to uncertainty.  It captures the following insight from "long-run risk"  models:  investor concerns of *long-term uncertainty* impacts *short-term asset valuation.*  In contrast to the "long-run risk" literature, our analysis opens the door to a different interpretation.  Instead of aversion to risk it reflects an aversion to the misspecification of to models or simplified perspectives on macroeconomic dynamics.  

As we noted in Remark {prf:ref}`remark:zero`,  $\left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right) W_{t+1}$ is  approximately the martingale component of the logarithm of consumption when  $\lambda$ is close to one.  In  {ref}`section:var` we showed that the variance of this component is challenging to estimate, a point originally made by {cite}`HansenHeatonLi:2008`.  This finding is part of the reason that we find it important to step back from rational expectations and limit investors confidence in the models they use for decision making.  

%Reference to {ref}`section:var`
````








(sec:second_order)=
### Order two


Differentiating equation {eq}`value_recur5` a second time gives:
```{math}
{\widehat V}_t^2  = 
(1 - \lambda)  {\widehat C}_t^2 +  \lambda  {\widehat R}_t^2  + (1 - \rho) (1 - \lambda) \lambda  \left( {\widehat R}_t^1 - {\widehat C}_t^1 \right)^2.
```
Equivalently,
```{math}
{\widehat V}_t^2 - {\widehat C}_t^2 = \lambda \left( {\widehat R}_t^2 - {\widehat C}_t^2 \right) +  (1 - \rho) (1 - \lambda) \lambda  \left( {\widehat R}_t^1 - {\widehat C}_t^1 \right)^2.
```

Rewrite transformations {eq}`tilde_transform` and {eq}`tilde_R` as 
```{math}
{\sf q} {\widetilde V}_t = {\widehat V}_t - {{\widehat V}_t^0} 
```
```{math}
{\sf q} {\widetilde R}_t =  {\widehat R}_t - {{\widehat V}_{t+1}^0}. 
```
Differentiating twice with respect to ${\sf q}$ and evaluated at ${\sf q} = 0$ gives:
```{math}
\left. 2 \frac{d} {d{\sf q}}{\widetilde V}_t  + {\sf q}  \frac{d^2} {d{\sf q}^2} {\widetilde V}_t \right\vert_{{\sf q} = 0} =  2{\widetilde V}_t^1   = {\widehat V}_t^2.
```
```{math}
\left. 2 \frac{d} {d{\sf q}}{\widetilde R}_t  + {\sf q}  \frac{d^2} {d{\sf q}^2} {\widetilde R}_t \right\vert_{{\sf q} = 0} =  2{\widetilde R}_t^1 =  {\widehat R}_t^2.
```
Differentiating {eq}`tilde_R` with respect to ${\sf q}$ gives:
```{math}
{\frac{d{\widetilde R}_t} {d {\sf q}}} =   \frac{{\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right] {\frac{d{\widetilde V}_{t+1}} {d {\sf q}}} \mid {\mathfrak A}_t \right)}{{\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right]  \mid {\mathfrak A}_t \right)}.
```
and thus
```{math}
:label: recur_update
{\widehat R}_t^2 = 2 {\widetilde R}_t^1 = 2 E \left( N_{t+1}^0  {\widetilde V}^1_{t+1} \mid {\mathfrak A}_t \right) =  E \left( N_{t+1}^0  {\widehat V}^2_{t+1} \mid {\mathfrak A}_t \right),
```
where subtracting ${\widehat C}_t^2$ from ${\widehat R}_t^2$ gives:
```{math}
:label: second_order_risk
{\widehat R}_t^2 - {\widehat C}_t^2 = E \left[ N_{t+1}^0 \left( {\widehat V}^2_{t+1} - {\widehat C}_{t+1}^2\right) + \left( {\widehat C}^2_{t+1} - {\widehat C}_{t}^2 \right)  \mid {\mathfrak A}_t \right] .
```
Substituting this formula into {eq}`recur_update` gives:
```{math}
:label: second_recursive_update 
\begin{align}
{\widehat V}_t^2 - {\widehat C}_t^2 =  \hspace{.2cm} & \lambda {\mathbb E} \left( N_{t+1}^0  \left[\left({\widehat V}^2_{t+1} - {\widehat C}^2_{t+1} \right) + \left( {\widehat C}^2_{t+1} -  {\widehat C}^2_{t}\right) \right] \mid {\mathfrak A}_t \right) \cr
\hspace{.2cm} & +  (1 - \rho) (1 - \lambda) \lambda  \left( {\widehat R}_t^1 - {\widehat C}_t^1 \right)^2.
\end{align}
```
Even if the second-order contribution to the consumption process is zero, there will be nontrivial adjustment to the approximation of ${\widehat V} - {\widehat C}$ because $\left( {\widehat R}^1 - {\widehat C}^1 \right)^2$ is different from zero. This term vanishes when $\rho = 1,$ and its sign will be different depending on whether $\rho$ is bigger or smaller than one.


## Stochastic discount factor approximation

We approximate $\left[{\widehat S}_{t+1} - {\widehat S}_{t} \right]$ in formula {eq}`eqn:sdf50` as
```{math}
{\widehat S}_{t+1} - {\widehat S}_{t}  \approx \left[{\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right] + \left[{\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right] + {\frac 1 2} 
\left[{\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right]
```
where
```{math}
{\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \eqdef  \log \beta - \rho \eta_c^0 
``` 
```{math}
{\widehat S}_{t+1}^1 - {\widehat S}_{t}^1 \eqdef - {\widehat C}_{t+1}^1 + {\widehat C}_{t}^1 + (\rho - 1) \left( {\widehat V}_{t+1}^1 - {\widehat C}^1_{t+1} \right)  - 
(\rho - 1) \left( {\widehat R}_t^1 - {\widehat C}_t^1  \right)
```
```{math}
{\widehat S}_{t+1}^2 - {\widehat S}_{t}^2 \eqdef - {\widehat C}_{t+1}^2 + {\widehat C}_{t}^2 + (\rho - 1) \left( {\widehat V}_{t+1}^2 - {\widehat C}^2_{t+1} \right)  - 
(\rho - 1) \left( {\widehat R}_t^2 - {\widehat C}_t^2  \right)
```

We now consider two different approaches to approximating $N_{t+1}^*$.  

###  Approach 1

Write
```{math}
N_{t+1}^* = \frac { \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right]}{{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right] \mid {\mathfrak A}_t \right]} 
```
```{math}
= \frac { \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right]}{ \exp\left[(1-\gamma_o) {\widetilde R}_t\right] } 
```
Form the ``first-order'' approximation:
```{math}
\log N_{t+1}^*  \approx (1 - \gamma_o) \left[\left( {\widetilde V}_{t+1}^0 - {\widetilde R}_t^0  \right) +   {\sf q} \left( {\widetilde V}_{t+1}^1 - {\widetilde R}_t^1 \right)\right] 
```
```{math}
= (1 - \gamma_o)\left[ \left( {\widehat V}_{t+1}^1 - {\widehat R}_t^1 \right) +  {\frac {\sf q}  2}  \left( {\widehat V}_{t+1}^2 - {\widehat R}_t^2 \right) \right]
```

We  combine a first-order approximation of $\log N_{t+1}^*$ with a second-order approximation of 
${\widehat S}_{t+1} - {\widehat S}_t$:
```{math}
\log S_{t+1} -  \log S_t \approx (1 - \gamma_o)\left[ \left( {\widehat V}_{t+1}^1 - {\widehat R}_t^1 \right) +  {\frac 1 2}  \left( {\widehat V}_{t+1}^2 - {\widehat R}_t^2 \right) \right] 
```
```{math}
+ \left({\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right) + \left({\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right)  + \frac 1 2 \left({\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right) 
```
which preserves the quadratic approximation of $\log S_{t+1} -  \log S_t.$. Note that If we were to use a second-order approximation of $N_{t+1}^*$, it would push us outside the class of exponentially quadratic stochastic discount factors. 


###  Approach 2



Next consider an alternative modification of Approach 1 whereby:
```{math}
:label: restricted_probability
\begin{align}
N_{t+1}^* \approx & \frac { \exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\widetilde V}_{t+1}^1 \right)\right]} 
{{\mathbb E}\left(\exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\widetilde V}_{t+1}^1 \right)\right] \mid {\mathfrak A}_t \right)} \cr
= & \frac { \exp\left[(1-\gamma_o)\left[ \left({\widehat V}_{t+1}^1 - {\widehat R}_t^1\right) + {\frac 1 2} \left( {\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right) \right]\right]} 
{{\mathbb E}\left(
\exp\left[(1-\gamma_o)\left[ \left({\widehat V}_{t+1}^1 - {\widehat R}_t^1\right) + {\frac 1 2} \left( {\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right) \right]\right]
\mid {\mathfrak A}_t \right)} \cr
\eqdef & {\widetilde N}_{t+1}
\end{align}
```
and $\log {\widetilde N}_t$ is used in conjunction with 
```{math}
\left({\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right) + \left({\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right)  + \frac 1 2 \left({\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right) . 
```
By design, this approximation of $N_{t+1}^*$ will have conditional expectation equal to one in contrast to the approximation used with Approach 1.   With a little bit of algebraic manipulation, it may be shown that this approximation induces a distributional change for $W_{t+1}$ with a conditional mean that is affine in $X_{t+1}$ and an altered conditional variance matrix that is constant over time.  

To understand better this choice of approximation, consider the family of random variables (indexed by ${\sf q}$) 
```{math}
(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\sf q} {\widetilde V}_{t+1}^1\right) 
- \log {\mathbb E} \left( \exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 +  {\sf q} {\widetilde V}_{t+1}^1\right) \right] \mid {\mathfrak A}_t \right) .
```
The corresponding family of exponentials has conditional expectation one and the ${\sf q} = 1$ member is the proposed approximation for $N_{t+1}^*.$ Differentiate the family with respect to ${\sf q}$:
```{math}
{\widetilde V}_{t+1}^1 - \frac {{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1}^0  \right] {\widetilde V}_{t+1}^1 \mid {\mathfrak A}_t \right)}
{{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1}^0  \right]  \mid {\mathfrak A}_t \right)} 
 =
{\widetilde V}_{t+1}^1 - {\widetilde R}_{t}^1.
```
Thus this family of random variables has the same first-order approximation in ${\sf q}$ as the one we derived previously for $\log N_{t+1}^*$.  

As a change of probability measure, this approximation will induce state dependence in the conditional mean and will alter the covariance matrix of the shock vector. We find this approach interesting because it links back directly to the outcome of the robustness formulation we described in [Section 3.1](sec:robust). 

## Solving a planner's problem with recursive utility

The {cite}`bansalyaron2004` example along with many others building connections between the macro economy and asset value take aggregate consumption as pre-specified. As we open the door to a richer collection of macroeconomic models, it becomes important to entertain more endogeneity, including investment and other variables familiar to macroeconomics.

Write a triangular system with stochastic growth as:
```{math}
:label: triangle
\begin{align}
X_{t+1} \left( \mathsf{q} \right) =  & \hspace{.2cm} \psi^x \left[D_t \left( \mathsf{q} \right), X_{t} \left( \mathsf{q} \right), {\sf q} W_{t+1}, {\sf q} \right] \cr
\log G_{t+1} \left( \mathsf{q} \right) - \log G_t \left( \mathsf{q} \right) = & \hspace{.2cm} \psi^g \left[ D_t \left( \mathsf{q} \right), X_{t} \left( \mathsf{q},
 \right), {\sf q} W_{t+1}, {\sf q}  \right],
 \end{align}
```
where $D_t$ is a date $t$ decision vector for the planner. Define ${\widehat G}_t = \log G_t.$ In addition, we impose 
```{math}
:label: output
{\widehat C}_t \left( \mathsf{q} \right) = \kappa \left[D_t \left( \mathsf{q} \right), X_{t} \left( \mathsf{q}  \right) \right] + {\widehat G}_t \left( \mathsf{q} \right).
```
In what follows for notational convenience, we will leave the ${\sf q}$ implicit.  
We use homogeneity to rewrite the utility recursion:
```{math}
{\widehat V}_t - {\widehat G}_t = \frac{1}{1-\rho} \log \left[ (1 - \beta) \exp\left[ (1-\rho) ({\widehat C}_t - {\widehat G}_t) \right] + \beta \exp\left[ (1-\rho) ({\widehat R}_t - {\widehat G}_t) \right] \right]  
```
where
```{math}
{\widehat R}_t - {\widehat G}_t = \frac{1}{1-\gamma} \log \left( {\mathbb E} \left[ \exp\left( (1 - \gamma) \left[ \left( {\widehat V}_{t+1} - {\widehat G}_{t+1} \right) + \left( {\widehat G}_{t+1} - {\widehat G}_t \right) \right] \right) \mid {\mathfrak A}_t \right] \right).
```
The approximation formulas, {eq}`first_order_risk`, {eq}`first_recursive_update`, {eq}`second_order_risk`, {eq}`second_recursive_update`,  that we deduced previously for ${\widehat V}_t - {\widehat C}_t$ and ${\widehat R}_t - {\widehat C}_t$ have and immediate counterparts for 
${\widehat V}_t - {\widehat G}_t$ and  ${\widehat R}_t - {\widehat G}_t$, from which can build an approximation of 
$ {\widehat R}_t - {\widehat V}_t.$ 



The first-order conditions for $D$ are:
```{math}
\begin{align*}
&(1-\beta) \exp\left[ (1- \rho) \left({\widehat C}_t - {\widehat G}_t \right) \right] \kappa_d( D_t, X_t) 
\cr &+ 
\beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
{\mathbb E} 
\left(
\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{d'}^x ( D_t, X_t, W_{t+1} )' MX_{t+1} \mid {\mathfrak A}_t \right) \cr
&+  \beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
{\mathbb E} 
\left(
\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat G}_t \right) \right] \psi_{d'}^g ( D_t, X_t, W_{t+1} ) \mid {\mathfrak A}_t \right)  \cr
& = 0.
\end{align*}
```
where $MX_{t+1}$ is the co-state, or the partial derivation of next period's value function evaluated at the next period's state vector. 
%Multiply both sides by $ \exp\left[ (1 - \rho) \left({\widehat V}_t - {\widehat G}_t \right) \right] $ to get
%```{math}
%\begin{align*}
%&(1-\beta) \exp\left[ (1- \rho) \left({\widehat C}_t - {\widehat G}_t \right) \right] \kappa_d( D_t, X_t)
%\cr
%& + 
%\beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
%{\mathbb E} 
%\left(
%\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{d'}^x ( D_t, X_t, W_{t+1} )' MX_{t+1} \mid {\mathfrak A}_t \right) \cr
%&+  \beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
%{\mathbb E} 
%\left(
%\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{d}^g ( D_t, X_t, W_{t+1} ) \mid {\mathfrak A}_t \right)  \cr
%& = 0.
%\end{align*}
%```
In addition, we solve a forward-looking co-state equation given by 
```{math}
\begin{align*}
& MX_t =  (1 - \beta) \exp\left[ (1 - \rho) \left({\widehat C}_t - {\widehat G}_t \right) \right] \kappa_{x}(D_t, X_t)   \cr
& + \beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
{\mathbb E} 
\left(
\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{x'}^x ( D_t, X_t, W_{t+1} )' MX_{t+1} \mid {\mathfrak A}_t \right) \cr
& +  \beta \exp\left[ (1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right) \right] 
{\mathbb E} 
\left(
\exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{x}^g ( D_t, X_t, W_{t+1} ) \mid {\mathfrak A}_t \right).  
\end{align*} 
```
Recall that $N_{t+1}^* =  \exp\left[ (1-\gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right].$  This allows us to compute the conditional expectations using the uncertainty adjusted transition probability.  

The solution method we use  finds  $MX_t$ and $D_t$ as time-invariant functions of $X_t$ that make the dynamic system behave in a stochastically stable manner.  

### An example economy with long-run uncertainty

Consider an AK model with recursive utility and adjustment costs.  

The exogenous state dynamic capture both long-run uncertainty in the mean growth rate and the overall volatility in the economy.  
```{math}
:label: equation2
\begin{align}
Z_{1,t+1} -
Z_{1,t}=   \hspace{.2cm} & - \epsilon \nu_{1} Z_{1,t}  +  \sqrt{\epsilon}\exp\left(\frac 1 2  {Z}_{2,t}\right) \sigma_1  W_{t+1} \cr 
{Z}_{2,t+1} - {Z}_{2, t} =  \hspace{.2cm} &   - \epsilon \nu_2 \left[ 1 - {\mu_2} \exp\left( - {Z}_{2,t} \right) \right]  \cr \hspace{.2cm}&  - {\frac \epsilon 2} |\sigma_2|^2  \exp\left( - {Z}_{2,t} \right)   
+ \sqrt{\epsilon} \exp\left( - {\frac 1 2}  {Z}_{2,t} \right) \sigma_2 W_{t+1}. \end{align} 
```
We include the $\epsilon$  to allow for a small time interval approximation of the continuous-time system.  A small value of $\epsilon$ is associated a time scaling in which a unit time interval is small, say months instead of years.      The state variable, $Z_{2,t}$ is included to capture stochastic volatility.  The discrete-time dynamics for $\{\exp(Z_{2,t}) \}$ approximate a continuous-time version of what is called a square root process  due to Feller.  Let
```{math}
X_t = \begin{bmatrix} Z_{1,t} \cr Z_{2,t} \end{bmatrix}.
```

With these exogenous dynamics, we obtain the following zero and first-order approximations: 
```{math} 
Z_{1,t}^0 = 0 \hspace{.4cm} \exp\left(Z_{2,t}^0\right)  =  \mu_2,
```
and 
```{math}
\begin{align}
Z_{1,t+1}^1 - Z_{1,t}^1 = & \hspace{.2cm}  - \epsilon \nu_1 Z_{1,t} + \sqrt{\epsilon \mu_2} \sigma_1 W_{t+1} \cr 
Z_{2,t+1}^1 - Z_{2,t}^1 = & \hspace{.2cm} - \epsilon \nu_2 Z_{2,t} + \sqrt{\frac \epsilon {\mu_2}} \sigma_2 W_{t+1} . 
\end{align}
```


We impose the resource constraint:
```{math}
C_t + I_t = \alpha K_t .
```
The endogenous state dynamics are given by:
```{math}
%:label: equation1
\begin{align} 
{{\widehat K}_{t+1}} - {\widehat K}_t = \hspace{.2cm} & \epsilon \left[ {\frac 1 \zeta}  \log \left( 1 + \zeta \frac{I_t}{K_t} \right) + \nu_k Z_{1,t} - \iota_k \right] \cr \hspace{.2cm} &  - \frac{\epsilon}{2}|\sigma_k|^2 \exp\left( {Z}_{2,t} \right)  
 +  \sqrt{\epsilon} \exp\left(\frac 1 2 { Z}_{2,t}  \right) {\sigma}_k  W_{t+1} 
  \end{align} 
```
where ${\widehat K}_t = \log K_t = {\widehat G}_t$.  The planner choice variable  $D_t =   \frac{I_t}{K_t}$.  
Rewrite the current-period  resource constraint as:  
```{math}
:label: equation3
{\widehat C}_t - {\widehat G}_t = \log \left( \alpha - D_t \right) .
```
The first-order conditions for the investment-capital ratio imply that  
```{math}
:label: equation4
\begin{align} 
0 = &\hspace{.2cm} - \exp\left[(1 - \rho) \left({\widehat C}_t - {\widehat G}_t \right)\right] \left(\frac {1-\beta}{\alpha - D_t}\right) 
\cr &\hspace{.2cm}  +  \exp\left[(1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right)\right] \left(\frac{\beta \epsilon }{1 + \zeta D_t} \right)
\end{align}    
```
To support a continuous-time approximation,  suppose $\beta = \exp(-\epsilon \delta)$ where $\delta$ is a the subjective rate of discount.   Notice that these first-order conditions do not depend on the co-state process $\{MX_t : t \ge 0 \}$.  We may solve this planner's problem using {eq}`equation4` and up-dating the continuation value processes and its uncertainty-adjusted counterpart until convergence.  When $\rho = 1$, ${\widehat R}_t - {\widehat G}_t$ and ${\widehat C}_t - {\widehat G}_t$ drop out of the first-order conditions and the solution for $D_t$ is constant.  When  $\rho \ne 1$, $D_t$ depends on the exogenous state $X_t$.  


### First-order approximation when $\rho=1$

%The Markov process governing the predictable component of macroeconomic growth is scalar in the {cite}`bansalyaron2004` analysis. Motivated by empirical evidence, {cite}`HansenHeatonLi:2008` study an extension of this model where $X^1$ is a vector autoregression. 

We obtain the following zero and first-order approximations for the  exogenous dynamics:
```{math} 
Z_{1,t}^0 = 0 \hspace{.4cm} \exp\left(Z_{2,t}^0\right)  =  \mu_2,
```
and 
```{math}
\begin{align}
Z_{1,t+1}^1 - Z_{1,t}^1 = & \hspace{.2cm}  - \epsilon \nu_1 Z_{1,t}^1 + \sqrt{\epsilon \mu_2} \sigma_1 W_{t+1} \cr 
Z_{2,t+1}^1 - Z_{2,t}^1 = & \hspace{.2cm} - \epsilon \nu_2 Z_{2,t}^1 + \sqrt{\frac \epsilon {\mu_2}} \sigma_2 W_{t+1} . 
\end{align}
```

The first-order conditions for $D$ are:
```{math}
\frac {1-\beta}{\alpha -D} + \frac {\beta \epsilon}{1 + \zeta D} = 0.
```
Solving for $D$ gives:
```{math} 
D^* = \frac {(\beta - 1)  + \beta \epsilon \alpha} {\beta \epsilon+ (1-\beta) \zeta } ,
```
which is independent of the state, as should be expected since $\rho = 1.$ From the capital evolution it follows from the order zero approximation is
```{math}
{\widehat K}_{t+1}^0 - {\widehat K}_t^0 = \epsilon \left[ \frac 1 \zeta \log  \left( 1 + \zeta D^*\right) - \iota_k \right]. 
```
The order one approximation is then:
```{math}
{\widehat K}_{t+1}^1 - {\widehat K}_t^1  =  \epsilon \nu_k Z_t^1 + \sqrt{\epsilon \mu_2} \sigma_k W_{t+1} . 
```
Stochastic volatility, as in the  {cite}`bansalyaron2004` model of consumption dynamics,  will be present in the second-order approximation.  

### Second-order approximation when $\rho = 1$

We next consider the second-order approximations. The second-oder approximation for $\{ Z_{2,t}\}$ does not contribute to the planner's solution or to the implied shadow prices and thus we drop it from the analysis.  For the remaining two state variables, we find that 

\begin{align*}
Z_{1,t+1}^2 - Z_{1,t}^2 = \hspace{.2cm} &   - \epsilon \nu_1 Z_{1,t}^2 - \epsilon \mu_2 |\sigma_1|^2   + \sqrt{\epsilon \mu_2} Z_{2,t}^1 \sigma_1 W_{t+1}  \cr
{\widehat K}_{t+1}^2 - {\widehat K}_t^2 = \hspace{.2cm} & -  \epsilon \mu_2 |\sigma_k|^2  + \sqrt{\epsilon \mu_2} Z_{2,t}^1 \sigma_k W_{t+1}
\end{align*} 
where we previously noted that $\{ Z_{2,t}^1\}$ evolves as first-order autoregression.  

The combined approximation for ${\sf q}=1$ uses:
\begin{align*}
Z_{1,t+1} - Z_{1,t} \approx & \hspace{.2cm}  - \epsilon \nu_1 Z_{1,t}  
\cr & \hspace{.2cm}  - \frac {\epsilon \mu_2 |\sigma_1|^2}  2 + \sqrt{\epsilon \mu_2} W_{t+1} + \frac {\sqrt{\epsilon \mu_2}    Z_{2,t} } 2 \sigma_1 W_{t+1} \cr
Z_{2,t+1} - Z_{2,t} \approx  & \hspace{.2cm} - \epsilon \nu_2 Z_{2,t} + \sqrt{\frac \epsilon {\mu_2}} \sigma_2 W_{t+1} \cr
{\widehat K}_{t+1} - {\widehat K}_t \approx  & \hspace{.2cm}  \epsilon \left[ \frac 1 \zeta \log  \left( 1 + \zeta D^*\right) - \iota_k \right]
+ \epsilon \nu_k Z_{1,t}  \cr & \hspace{.2cm}  - \frac {\epsilon \mu_2 |\sigma_k|^2}  2  + \sqrt{\epsilon \mu_2} \sigma_kW_{t+1}     + \frac {\sqrt{\epsilon \mu_2} Z_{2,t}}2 \sigma_k W_{t+1}.
\end{align*}

The approximate dynamics for the exogenous states remains the same for $\rho \ne 1,$ but the solution for $D^*$ becomes state dependent and the approximate dynamic evolution for capital is altered.  








%
%### Shock elasticities
%
%We use the shock elasticities to explore pricing implications of this recursive utility specification. We conduct this exploration using the original parameter calibration in {cite}`bansalyaron2004`. These elasticities and their relation to impulse-response functions introduced first to macroeconomics by {cite}`frisch:1933` is described in {cite}`BorovickaHansenScheinkman:2014`. In what follows, we use exponential/linear/quadratic implementation by {cite}`borovickahansen14` and by {cite}`borovivcka2016term`. 
%
%{numref}`fig:expo` gives the shock exposure elasticities for consumption to each of the three shocks. This can be interpreted as nonlinear local impulse responses for consumption (in levels not logarithms). The elasticities for the growth rate shock and the stochastic volatility shock start small and increase over the time horizon as dictated by the persistence of the two exogenous state variable processes. The elasticities for the direct shock to consumption are flat over the horizon as to be expected since the shock directly impacts log consumption in a manner that is permanent. Notice that while elasticities for the volatility shock are different from zero, their contribution is much smaller than the other shocks.
%
%```{figure} ../images/recursive/Expo.png
%:name: fig:expo
%:width: 1000px
%Exposure elasticities for three shocks. The time scale is in months. 
%```
%
%
%```{prf:footnote}
% :label: footnote2
%It is notable that we are looking at levels and not logarithms of consumption. The local impulse response for the logarithms of consumption is in fact zero for the stochastic volatility shock.
%```
%
%Nevertheless, for this {cite}`bansalyaron2004` calibration of the long run risk model, stochastic volatility induces state dependence in the elasticities for growth rate and consumption shocks as reflected by quantiles given in the figures.  
%
%
%{numref}`fig:price1` gives the corresponding shock price elasticities for $\rho = 2/3$ and $\gamma = 10$. The recursive utility preferences are forward looking as reflected by the continuation-value contribution to the one-period increment to the stochastic discount factor process as given in {eq}`eqn:sdf50`. This forward-looking contribution is reflected in shock price elasticities that are now flat for both the growth rate shock and the shock to stochastic volatility. The magnitudes are substantially higher for the shock-price elasticities. While the relative magnitudes are very different, the shock price elasticities are much smaller than the other elasticities.  
%
%```{prf:footnote}
% :label: footnote3
%We normalized the stochastic volatility shock $\sigma_x^2$ to be negative implying that a positive shock reduces the stochastic volatility state variable. Under this normalization, the shock price elasticities are positive.
%```
%
%```{figure} ../images/recursive/Price_006.png
%:name: fig:price1
%:width: 1000px
%Price elasticities for three shocks. $\rho = 2/3, \gamma = 10,$ $\beta = .998.$  The time scale is in months.
%```
%
%{numref}`fig:price2` and {numref}`fig:price3` provide the analogous plots for $\rho = 1, 1.5.$ The shock price elasticities are very similar given these modes increases in $\rho$. It is evidently the risk aversion parameter $\gamma = 10$ that is important for determining the magnitude of these elasticities. {numref}`fig:price4` sets $\rho=\gamma = 10$ which corresponds to preferences that are time separable. The forward-looking component to the stochastic discount factor is shut down as is evident from formula {eq}`eqn:sdf50`. Now the shock price elasticities and shock exposure elasticities show a very similar trajectory except that the shock price elasticities are about ten times larger. The stochastic volatility shock price is increased by about seventy-five times.
%
%```{figure} ../images/recursive/Price_010.png
%:name: fig:price2
%:width: 1000px
%Price elasticities for three shocks. $\rho = 1, \gamma = 10,$ $\beta = .998.$  The time scale is in  months.  
%```
%
%```{figure} ../images/recursive/Price_015.png
%:name: fig:price3
%:width: 1000px
%Price elasticities for three shocks. $\rho = 3/2, \gamma = 10,$ $\beta = .998.$  The time scale is in  months. 
%```
%
%```{figure} ../images/recursive/Price_100.png
%:name: fig:price4
%:width: 1000px
%Price elasticities for three shocks. $\rho = 10, \gamma = 10,$ $\beta = .998.$ The time scale is in  months.  
%```
%

## Solving models

 In this section, we briefly describe one way to extend the approach that builds directly on previous second-order approaches of {cite}`KimKimSchaumbergSims:2008`, {cite}`schmitt2004solving`, and {cite}`LombardoUhlig:2018`. While such methods should not be viewed as being generically applicable to nonlinear stochastic equilibrium models, we find them useful pedagogically and often as at least initial steps to understanding models that are arguably “smooth.” See {cite}`PohlSchmeddersWilms:2018` for a careful study of nonlinearity in asset pricing models with recursive utility.[^pohl]

We implement these methods for second-order approximation using the following steps.

1. Solve for ${\sf q}=0$ deterministic model.

2.  Take as given first and second-order approximate solutions for ${\widehat C}_t - {\widehat G}_t$ and ${\widehat G}_{t+1} - {\widehat G}_t.$. Solve for the approximate solutions for ${\widehat V}_t - {\widehat G}_t,$ ${\widehat V}_{t+1} - {\widehat R}_t$ and $N_{t+1}.$ 


3. Compute the first-order expansion and solve the resulting equations following the previous literature for $D_t,$ ${\widehat C}_t - {\widehat G}_t,$ and ${\widehat G}_{t+1} - {\widehat G}_t.$ When constructing these equations, use expectations computed using the probabilities induced by $N_{t+1}^0$.  Substitute the first-order approximation for ${\widehat R}_t - {\widehat G}_t.$


4. Compute the second-order expansion and solve the resulting equations following the previous literature. Again use the expectations induced by $N_{t+1}^0$. In addition, make another recursive utility adjustment expressed in terms the approximations of  ${\widehat R}_t - {\widehat G}_t.$

5.  Return to step 2, and repeat until convergence.  


Initialize this algorithm by solving the $\gamma_o=1$ and $\rho = 1,$ which can be solved without iteration.  


 

See the Appendices that follow for more details and formulas to use in the solution method.

[^pohl]: {cite}`PohlSchmeddersWilms:2018` provide examples of when log-linear or local methods of computation fail to provide good approximations.



As a second approach we iterate over an $N_{t+1}^*$ given by formula {eq}`restricted_probability` approximation restricted to induce an alternative probability distribution Call the approximation ${\widetilde N}_{t+1}$ with an induced distribution for $W_{t+1}$ that is normal with conditional mean ${\tilde \mu}_t$ and covariance matrix ${\widetilde \Sigma}$. 





While we discussed the approximation for resource allocation problems with recursive utility, there is a direct extension of this approach to solve a  general class stochastic equilibrium models by stacking a system of expectational-type equations expressed in part using the recursive utility stochastic discount factor that we derived.  For resource allocation problems, we expressed the first-order conditions for the planner in utility units, which simplified some formulas. Equilibrium models not derived from a planner's problems typically use  stochastic discount factors expressed in consumption units when representing investment choices.  The approximation methods described in this chapter have a direct extension to such models.
































---

## Appendix A: Solving the planner's problem


Consider the equation:
```{math}
Q_t {\mathbb E} \left( N_{t+1} H_{t+1} \mid {\mathfrak A}_t \right)  + L_{t}   = 0 
```
where
\begin{align*}
Q_t \eqdef & \hspace{.2cm}  \beta \exp\left[(1-\rho) \left( {\widehat R}_t - {\widehat G}_t \right) \right] \cr \cr 
H_{t+1}  \eqdef & \hspace{.2cm} \begin{bmatrix} \psi_{d'}^x (D_t, X_t, W_{t+1})'MX_{t+1} +  \psi_{d'}^g (D_t, X_t, W_{t+1}) \cr \psi_{x'}^x (D_t, X_t, W_{t+1})'MX_{t+1} +  \psi_{x}^g (D_t, X_t, W_{t+1}) \cr
\end{bmatrix} \cr \cr
L_t \eqdef & \hspace{.2cm} \begin{bmatrix} (1-\beta) \exp \left[ (1-\rho) \left({\widehat C}_t - {\widehat G}_t \right) \right] \kappa_d(D_t,X_t)   \cr 
(1-\beta) \exp \left[ (1-\rho) \left({\widehat C}_t - {\widehat G}_t \right) \right] \kappa_x(D_t, X_t) - MX_t \end{bmatrix} 
.
\end{align*}

Solve for $ {\widehat C}_t - {\widehat G}_t, D_t, MX_t$ as a function of $X_t$.    We include the state dynamics {eq}`triangle` and the output constraint {eq}`output` when computing a solution.  The objects: ${\widehat C} - {\widehat G}, D, MX$ are sometimes referred to as jump variables since  we not impose initial conditions as part of a solution, in contrast to $X$.  

Our solution will entail an iteration.  We will impose a specification for $Q$ and $N$ and find an approximate solution for the dynamical system.  Then given this solution, we will compute a new implied solution for $Q$ and $N$.  We then iterate this until we achieve numerical convergence.   We use second-order approximations for both steps.  


###  $Q$ derivatives

To construct a candidate $Q$, we use the following strategy.  Compute ${\widehat V}_t^0 - {\widehat G}^0_t$ and ${\widehat G}_{t+1}^0-{\widehat G}_t^0$ as part of the steady state, and form:
```{math}
Q^0_t \eqdef \beta \exp[(1-\rho) \left[\left( {\widehat V}_{t+1}^0 - {\widehat G}_{t+1}^0\right)  + \left( {\widehat G}_{t+1}^0 - {\widehat G}_t^0 \right) \right] 
```
Construct
```{math} 
{\widehat R}_t^1 - {\widehat G}_t^1 = \frac 1 {1-\gamma_o} {\mathbb E} \left( \exp \left[ (1 - \gamma_o) \left( {\widehat V}_{t+1}^1 - {\widehat G}_{t+1}^1\right) +
\left( {\widehat G}_{t+1}^1 - {\widehat G}_{t}^1\right) \right] \mid {\mathfrak A}_t \right),
```
and 
\begin{align*}
Q_t^1 \eqdef & \hspace{.2cm} \beta (1-\rho) Q_t^0 \left( {\widehat R}_t^1 - {\widehat G}_t^1\right) \cr
Q_t^2 \eqdef  & \hspace{.2cm} \beta (1-\rho)^2 Q_t^0 \left( {\widehat R}_t^1 - {\widehat G}_t^1\right)^2 + 
\beta (1-\rho) Q_t^0 \left( {\widehat R}_t^2 - {\widehat G}_t^2\right).  
\end{align*} 

We rewrite {eq}`first_recursive_update` as:
\begin{align*}
& {\widehat V}_t^1 - {\widehat G}_t^1 =  (1-\lambda) \left({\widehat C}_t^1 - {\widehat G}_t^1 \right)  \cr
& + \left( {\frac \lambda {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o) \left[  \left( {\widehat V}_{t+1}^1  - {\widehat G}_{t+1}^1 \right) + \left( {\widehat G}_{t+1}^1 - {\widehat G}_t^1\right) \right] \right] \mid {\mathfrak A}_t \right)
\end{align*}
and solve this equation forward by first computing  the $\gamma_o=1$ answer and then adjusting this answer  for $\gamma_o > 1$ analogous to the approach described in {prf:ref}`remark_solve_first`.   

Given this solution, form
```{math}
 {\widehat R}_t^1 - {\widehat G}_t^1 =  \left( \frac 1 \lambda \right) \left( {\widehat V}_t^1 - {\widehat G}_t^1\right)
- \left( \frac  {1-\lambda} \lambda \right)  \left({\widehat C}_t^1 - {\widehat G}_t^1 \right).
```

For the second-order, it follows from {eq}`second_recursive_update` that 
```{math}
\begin{align}
{\widehat V}_t^2 - {\widehat G}_t^2 =  \hspace{.2cm} & 
\lambda {\mathbb E} \left( N_{t+1}^0  \left[\left({\widehat V}^2_{t+1} - {\widehat G}^2_{t+1} \right) + \left( {\widehat G}^2_{t+1} -  {\widehat G}^2_{t}\right) \right] \mid {\mathfrak A}_t \right) \cr
\hspace{.2cm} & + (1- \lambda) \left({\widehat C}_t^2 - {\widehat G}_t^2 \right)  +  (1 - \rho) (1 - \lambda) \lambda  \left( {\widehat R}_t^1 - {\widehat G}_t^1 + {\widehat G}_t^1  - {\widehat C}_t \right)^2,
\end{align}
```
which we solve this equation forward under the $N^0$ implied change in probability measure.  Form:
```{math}
\begin{align}
{\widehat R}_t^2 - {\widehat G}_t^2 = & \left( \frac 1 \lambda \right)  \left({\widehat V}_t^2 - {\widehat G}_t^2\right) 
- \left( \frac {1 - \lambda}{\lambda} \right) \left({\widehat C}_t^2 - {\widehat G}_t^2 \right) \cr 
& -  (1 - \rho) (1 - \lambda)   \left( {\widehat R}_t^1 - {\widehat G}_t^1 + {\widehat G}_t^1  - {\widehat C}_t \right)^2.
\end{align} 
```

### $N$ derivatives

```{math}
N_{t+1}^0 \eqdef \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1}^0 - {\widetilde R}_t^0 \right) \right]   = \exp\left[  (1-\gamma_o) \left({\widehat V}_{t+1}^1 - {\widehat  R}_t^1\right) \right]
```
```{math}
\begin{align}
N_{t+1}^1 \eqdef  \hspace{.2cm} & \left.  \frac d {d{\sf q}}  \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1} - {\widetilde R}_t \right) \right] \right\vert_{{\sf q} =0} \cr  =  \hspace{.2cm} & (1 - \gamma_o)  N_{t+1}^0\left( 
{\widetilde V}_{t+1}^1  - {\widetilde R}_t^1 \right) \cr =  \hspace{.2cm} & \left(\frac {1-\gamma_o} 2 \right) N_{t+1}^0 \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right).
\end{align}
```
Form:
\begin{align*} 
 {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 = & \left({\widehat V}_{t+1}^1 -  {\widehat G}_{t+1}^1\right) + 
 \left({\widehat G}_{t+1}^1 -  {\widehat G}_{t}^1\right) - \left({\widehat R}_{t}^1 -  {\widehat G}_{t}^1\right) \cr
 {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 = & \left({\widehat V}_{t+1}^2 -  {\widehat G}_{t+1}^2\right) + 
 \left({\widehat G}_{t+1}^2 -  {\widehat G}_{t}^2\right) - \left({\widehat R}_{t}^2 -  {\widehat G}_{t}^2\right).
\end{align*} 
It may be directly verified that $N_{t+1}^1$ and $ N_{t+1}^2$  have conditional expectations equal to zero.  Express
```{math} 
:label: V-R2
\begin{align}
 {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 = &  \left({\frac 1 {1-\gamma_o}} \right) \left[\mu^0 \cdot ( W_{t+1} - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \cr
 {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2 = &  {\frac 1 2} \left(W_{t+1} - \mu^0 \right)'\Upsilon_2^2 \left(W_{t+1} - \mu^0 \right) -
 {\frac 1 2} \rm{tr}\left( \Upsilon_2^2\right) \cr & -   \left(W_{t+1} - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right).
 \end{align}
```
In producing these representations, we use that have conditional ${\widehat V}_{t+1}^2 - {\widehat R}_{t}^2$ and  $ {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2$ have mean zero under the conditional probability distribution induced by $N_{t+1}^0.$ 





## Appendix B: Approximation formulas  (approach one)

Consider the equation:
```{math}
Q_t {\mathbb E} \left( N_{t+1}  H_{t+1} \mid {\mathfrak A}_t \right)  + L_{t}   = 0 .
```
where $ \beta \exp \left[- \rho \left( \log C_{t+1} - \log C_t \right) \right] $ is absorbed into the construction of $ H_{t+1} $. This is the subsystem of the equations not including the state evolution equations.  

(sec2)= 
### Order zero 

The order zero approximation of the product: $ N_{t+1} Q_{t+1} H_{t+1} +  L_{t} $ is:
```{math}
Q_t^0 N_{t+1}^0  H_{t+1}^0 + L_t^0 = 0 
```
Thus the order zero approximate equation is:
```{math}
Q_t^0 {\mathbb E} \left[N_{t+1}^0 \left( H_{t+1}^0  \right)  \mid {\mathfrak A}_t \right] + L_{t+1}^0=  Q_t^0 H_{t+1}^0 + L_{t}^0 = 0
```
since $ N_{t+1}^0 $ has conditional expectation equal to one.  We add to this subsystem the  ${\sf q} = 0$ state dynamic equation inclusive of jump variables,  and we compute a stable steady state solution.   

(sec3)= 
### Order one 

The order one approximation of the product: $ Q_t N_{t+1}  H_{t+1}+  L_{t} $ is:
```{math}
Q_t^1 N_{t+1}^0   H_{t+1}^0  + Q_t^0 N_{t+1}^1 H_{t+1}^0 + Q_t^0 N_{t+1}^0  H_{t+1}^1 +  L_{t}^1 .
```
Thus the order one approximate equation is:
```{math}
\begin{align}
& Q_t^0 {\mathbb E} \left( N_{t+1}^1   H_{t+1}^0    + N_{t+1}^0 H_{t+1}^1  \mid {\mathfrak A}_t \right)
+ Q_t^1 H_{t+1}^0    + L_t^1 \\
& = Q_t^0 {\mathbb E} \left(  N_{t+1}^0 H_{t+1}^1  \mid {\mathfrak A}_t \right) + Q_{t}^1 H_{t+1}^0 +  L_t^1 \\
& = 0
\end{align}
```
where we used the implication that $ {\mathbb E} \left( N_{t+1}^1  \mid {\mathfrak A}_t \right) = 0.$  


### Order two 

The order two approximation of the product: $ N_{t+1} Q_{t+1} H_{t+1} + L_{t+1} $ is:
```{math}
\begin{align}
& Q_t^0 N_{t+1}^0 H_{t+1}^2 + Q_t^0 N_{t+1}^2H_{t+1}^0  + 
2 Q_t^0 N_{t+1}^1 H_{t+1}^1 + L_t^2 \\
+ 
& 2 N_{t+1}^1  Q_{t}^1 H_{t+1}^0  +   2 N_{t+1}^0 Q_{t}^1 H_{t+1}^1 + N_{t+1}^0 Q_{t}^2 H_{t+1}^0 
\end{align}
```
The terms $ Q_t^0 N_{t+1}^2H_{t+1}^0$ and $2 N_{t+1}^1  Q_{t}^1 H_{t+1}^0$ have conditional expectation equal to zero.   Thus the approximating equation is:
```{math}
\begin{align}
& Q_t^0 {\mathbb E} \left( N_{t+1}^0 H_{t+1}^2 \mid {\mathfrak A}_t \right) +
L_t^2 \cr & + 2Q_t^0 {\mathbb E}\left( N_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) + 
2 Q_t^1 {\mathbb E} \left(N_{t+1}^0  H_{t+1}^1 \mid {\mathfrak A}_t \right) + H_{t+1}^0 Q_t^2  {\mathbb E} \left(  N_{t+1}^0 \mid {\mathfrak A}_t \right) \cr
& = 0.
\end{align}   
```

To elaborate on  the contributions in the second line, 
express $ H_{t+1}^1 $  as 
```{math}
:label: Hplus1
H_{t+1}^1  = \Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right). 
```
Then
```{math}
\begin{align} 
2Q_t^0 {\mathbb E}\left( N_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = (1-\gamma)  
  Q_t^0 {\mathbb E} \left[ N_{t+1}^0   \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right)H_{t+1}^1 \mid{\mathfrak A}_t \right], \\ 
& = 2(\gamma_o - 1) Q_t^o \Theta_2^1 \left[  \Upsilon_1^2 X_t^1  + \Upsilon_0^2\right] \\
2 Q_t^1 {\mathbb E} \left(N_{t+1}^0  H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 
2 Q_t^0 \beta (1- \rho) \left({\widehat R}_t^1 - {\widehat G}_t^1 \right)\left( \Theta_0^1 + \Theta_1^1 X_t^1 \right), \\
H_{t+1}^0 Q_t^2 {\mathbb E} \left(  N_{t+1}^0   \mid {\mathfrak A}_t \right) & =  H_{t+1}^0 Q_t^0  \left[\beta (1-\rho)^2  \left( {\widehat R}_t^1 - {\widehat G}_t^1\right)^2 + 
\beta (1-\rho)  \left( {\widehat R}_t^2 - {\widehat G}_t^2\right)\right].
\end{align} 
```
The formula for the first of these terms follows from {eq}`V-R2` and {eq}`Hplus1`, along with fact the third central moments of normals are zero.    

%Denote the sum of the four terms in the second affine as ${\overline H}_t^2$. This random variable will be affine in 
%$ X_t^1 $, with a dynamic evolution determined by solving the first-order approximation.  Thus we write the subsystem of equations to be solved as:
%```{math}
%{\mathbb E} \left( N_{t+1}^0 H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\overline H}_t^2  = 0.
%```
We add to this second-order subsystem, the second-order approximation of the state dynamics inclusive of the jump variables.  We substitute in the solution for the first-order approximation for the jump variables into both the first and second-order approximate state dynamics.  In solving the second-order jump variable adjustment we use expectations induced by $ N_{t+1}^0 $ zero throughout under which $W_{t+1} $ is conditionally normally distributed with mean $ \mu^0 $ and covariance $ I $.

## Appendix C: Approximation formulas  (approach two)

In this approach we use the same order zero approximation.  For the order one approximation, we use the formula {eq}`restricted_probability` for ${\widetilde N}_{t+1},$ which approximates  $N_{t+1}^*,$ in conjunction with:
\begin{align}
Q_t^0 {\mathbb E} \left(  {\widetilde N}_{t+1} H_{t+1}^1  \mid {\mathfrak A}_t \right) + Q_{t}^1 H_{t+1}^0 +  L_t^1 = 0.
\end{align}
From formula {eq}`V-R2`, it follows that under the ${\widetilde N}_{t+1},$ induced change in probability, $W_{t+1}$
is normally distributed with conditional mean 
```{math}
\mu^0 +  (1-\gamma_o) \left(\Upsilon_2^2\right)^{-1} \left(\Upsilon_1^2X_t^1 + \Upsilon_0^2\right) ,
```
and conditional precision:
```{math}
(\gamma_o - 1) \Upsilon_2^2 + {\mathbb I}.  
```
For the order two approximation, we use:
```{math} 
Q_t^0 {\mathbb E} \left(  {\widetilde N}_{t+1}  H_{t+1}^2  \mid {\mathfrak A}_t \right) + L_t^2 
+ 2 Q_t^1  {\mathbb E} \left(  {\widetilde N}_{t+1} H_{t+1}^1  \mid {\mathfrak A}_t \right) + H_{t+1}^0 Q_t^2 = 0 . 
```


%### Steps for implementation 
%
%
%
%We implement these methods for second-order approximation using the following steps.
%
%
%
%1. Solve $Q_t^0$, $H_{t+1}^0,$ and $L_{t+1}^0$ for order zero state and jump variables. The outcome will be state invariant.
%
%2. Take as given a $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in representations {eq}`VminusR2`.
%
%3. Compute the first-order contribution to approximation by following the previous literature with expectations computed using the probabilities induced by $N_{t+1}^0$, which imply that $W_{t+1}$ has mean $\mu^0$. Express the solution as in {eq}`H1`.
%
%4. Compute the second-order contribution to the approximation by following the previous literature, again with the expectations induced by $N_{t+1}^0$.
%
%5. Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in representations {eq}`VminusR2` and return to {stepi}. Repeat until convergence.
%
%
%
%## Appendix: Approximation formulas for expectation equations (approach two)
%
%For this solution, we iterate over $N_{t+1}^*$ approximation. Call the approximation ${\widetilde N}_{t+1}$ with an induced distribution for $W_{t+1}$ that is normal with conditional mean ${\tilde \mu}_t$ and covariance matrix ${\widetilde \Sigma}$. This distribution is used in both the first-order and second-order contributions to the approximation. The conditional mean for ${\tilde \mu}_t$ is affine in $X_t^1$. The following delineates the changes that need to be made.
%
%### First-order adjustment
%Compute:
%```{math}
%\begin{align*}
%{\mathbb E} \left( {\widetilde N}_{t+1}  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right) = \hspace{.2cm} &
%(\rho - 1) 
%{\mathbb E} \left[ {\widetilde N}_{t+1} \left({\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
%= \hspace{.2cm}&\left({\frac {\rho -1} {1-\gamma_o}}\right)  \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] H_{t+1}^0 \cr
% \eqdef \hspace{2cm} & {\widetilde  H}_t^1.
%\end{align*}
%```
%Then the equation to be solved is:
%```{math}
%{\mathbb E} \left({\widetilde N}_{t+1}   H_{t+1}^1 \mid {\mathfrak A}_t \right) + L_{t}^1 + {\widetilde H}^1_t  = 0.
%```
%
%### Second-order adjustment
%```{math}
%:label: second_affine_again
%\begin{align*}
% 2 {\mathbb E} \left({\widetilde N}_{t+1} Q_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) =  \hspace{.2cm} & 2(\rho - 1) {\mathbb E}\left[ {\widetilde N}_{t+1}\left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right) \right] \mid {\mathfrak A}_t \right] \cr
% = \hspace{.2cm}  & 2 \frac {(\rho - 1)}{(1-\gamma_o)}  \Theta_2^1 {\widetilde \Sigma} \mu^0 
%\cr  + &2 \frac {(\rho - 1)}{(1-\gamma_o)} \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1\left( {\tilde \mu}_t - \mu^0\right)  \right] \cr
%{\mathbb E} \left(  {\widetilde N}_{t+1} Q_{t+1}^2 H_{t+1}^0 \mid {\mathfrak A}_t \right) =  \hspace{.2cm}& (\rho - 1)^2{\mathbb E} \left[{\widetilde N}_{t+1} \left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2  \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr 
%&  +   (\rho - 1) {\mathbb E} \left[{\widetilde N}_{t+1}  \left( 
%{\widehat V}_{t+1}^2  - {\widehat R}^2 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
%= \hspace{.2cm}& \left(\frac {1-\rho}{1-\gamma_o}\right)^2\left[{\mu^0}'{\widetilde \Sigma} \mu^0 + \left(\mu^0 \cdot {\tilde \mu}_t - {\frac 1 2} |\mu^0|^2 \right)^2\right]
%H_{t+1}^0 \cr 
% &  + \frac {(\rho - 1)} 2 \left[ \rm{tr}\left( \Upsilon_2^2{\widetilde \Sigma} - \Upsilon_2^2 \right)  
%  +  \left({\tilde \mu}_t - \mu^0\right)' \Upsilon_{2}^2  \left({\tilde \mu}_t - \mu^0\right)\right] H_{t+1}^0 \cr
%&  +  (\rho - 1) \left({\tilde \mu}_t  - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right) H_{t+1}^0.
%\end{align*}
%```
%
%Denote the sum of the two terms in {eq}`second_affine_again` as ${\widetilde H}_t^2.$ Then the equation to be solved is 
%```{math}
%{\mathbb E} \left( {\widetilde N}_{t+1} H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\widetilde H}_t^2  = 0.
%```
%
%### Updated recursive utility adjustments
%
%Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2,  \Upsilon_2^2$ used in representations {eq}`VminusR2`. Compute a new version of
%```{math}
%{\widetilde  N}_{t+1}  =  
%\frac 
%{\exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right]} 
%{{\mathbb E}\left( \exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right] \mid {\mathfrak A}_t \right)},
%```
%and deduce the implied ${\widetilde \mu}_t$ and ${\widetilde \Sigma}$. The conditional mean ${\tilde \mu}_t$ satisfies:
%```{math}
%{\widetilde \Sigma}^{-1} {\tilde \mu}_t =  \mu^0 +  {\frac {(1-\gamma_o)} 2} \left(\Upsilon_0^2 + \Upsilon_1^2 X_t^1 - \Upsilon_2^2 \mu^0 \right)
%```
%where the formula for ${\widetilde \Sigma}$ is
%```{math}
%\widetilde {\Sigma} =  \left[{\mathbb I} - {\frac {(1 - \gamma_o)}  2} \Upsilon_2^2\right]^{-1}.
%```
%
%With these adjustments, we iterate to convergence.
%
%
%
%## Appendix: Approximation formulas for the stochastic discount factor
%
%***This should be moved to a Jupyter notebook***
%
%<!-- Test reference {ref}`section:var` -->
%
%For the purposes of this appendix, write:
%```{math}
%:label: equation5
%\frac {S_{t+1}}{S_t} = N_{t+1}^* Q_{t+1} \beta \exp \left[- \rho \left( \log C_{t+1} - \log C_t \right) \right]
%```
%
%where
%```{math}
%:label: equation6
%N_{t+1}^* = \exp\left[  (1-\gamma_o) \left({\widetilde V}_{t+1} - {\widetilde  R}_t\right) \right]
%```
%```{math}
%Q_{t+1}^* = \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right]
%```
%are terms that are contributed by recursive utility.  
%
%(sec:N_derivatives)=
%### $N_{t+1}^*$ derivatives
%
%```{math}
%:N_{t+1}^0:
%N_{t+1}^0 \eqdef \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1}^0 - {\widetilde R}_t^0 \right) \right]   = \exp\left[  (1-\gamma_o) \left({\widehat V}_{t+1}^1 - {\widehat  R}_t^1\right) \right]
%```
%```{math}
%\begin{align*}
%N_{t+1}^1 \eqdef  \hspace{.2cm} & \left.  \frac d {d{\sf q}}  \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1} - {\widetilde R}_t \right) \right] \right\vert_{{\sf q} =0} \cr  =  \hspace{.2cm} & (1 - \gamma_o)  N_{t+1}^0\left( 
%{\widetilde V}_{t+1}^1  - {\widetilde R}_t^1 \right) \cr =  \hspace{.2cm} & \left(\frac {1-\gamma_o} 2 \right) N_{t+1}^0 \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right).
%\end{align*}
%```
%It may be directly verified that $N_{t+1}^1$ has conditional expectation equal to zero.
%
%### $Q_{t+1}^*$ derivatives
%
%```{math}
%Q_{t+1}^0 \eqdef \exp\left[(\rho - 1) \left({\widehat V}_{t+1}^0 - {\widehat R}_t^0 \right) \right]   = 1 
%```
%```{math}
%Q_{t+1}^1 \eqdef \left.  \frac d {d{\sf q}}  \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \right\vert_{{\sf q} =0}  = (\rho - 1) \left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) 
%```
%```{math}
%Q_{t+1}^2  \eqdef \left.  \frac {d^2}  {d{\sf q}^2}   \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \right\vert_{{\sf q} =0}   = 
%(\rho - 1)^2 \left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2 +  (\rho - 1) \left( 
%{\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right)
%```
%
%Express
%```{math} 
%:widehat{V}_{t+1}^1 - \widehat{R}_{t}^1:
% {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 =  {\frac 1 {1-\gamma_o}} \left[\mu^0 \cdot ( W_{t+1} - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] 
%```
%```{math} 
%:label: VminusR2
% {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2 = {\frac 1 2} \left(W_{t+1} - \mu^0 \right)'\Upsilon_2^2 \left(W_{t+1} - \mu^0 \right) - 
% {\frac 1 2} \rm{tr}\left( \Upsilon_2^2\right) +  \left(W_{t+1} - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right),
%```
%Recall from {eq}`recur_update` that $ {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2$ has mean zero under the probability distribution induced by $N_{t+1}^0,$ which is consistent with its representation in {eq}`VminusR2`.
%
%
%
%
%## Appendix: Approximation formulas for expectation equations (approach one)
%
%Consider the equation:
%```{math}
%{\mathbb E} \left( N_{t+1} Q_{t+1} H_{t+1} \mid {\mathfrak A}_t \right)  + L_{t}   = 0 .
%```
%where $ \beta \exp \left[- \rho \left( \log C_{t+1} - \log C_t \right) \right] $ is absorbed into the construction of $ H_{t+1} $. This is the subsystem of the equations not including the state evolution equations.  
%
%(sec2)= 
%### Order zero 
%
%The order zero approximation of the product: $ N_{t+1} Q_{t+1} H_{t+1} +  L_{t} $ is:
%```{math}
%N_{t+1}^0  H_{t+1}^0 + L_t^0 = 0 
%```
%where we have substituted $ Q_{t+1}^0 = 1 $.  
%Thus the order zero approximate equation is:
%```{math}
%{\mathbb E} \left[N_{t+1}^0 \left( H_{t+1}^0 + L_{t+1}^0 \right)  \mid {\mathfrak A}_t \right] =  H_{t+1}^0 + L_{t}^0 = 0
%```
%since $ N_{t+1}^0 $ has conditional expectation equal to one.  We add to this subsystem the  ${\sf q} = 0$ state dynamic equation inclusive of jump variables,  and we compute a stable steady state solution.   
%
%(sec3)= 
%### Order one 
%
%The order one approximation of the product: $ N_{t+1} Q_{t+1} H_{t+1}+  L_{t} $ is:
%```{math}
%N_{t+1}^1   H_{t+1}^0  + N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0 + N_{t+1}^0  H_{t+1}^1 +  L_{t}^1 
%```
%where we have substituted $ Q_{t+1}^0 = 1 $.  Thus the order one approximate equation is:
%```{math}
%\begin{align*}
%& {\mathbb E} \left( N_{t+1}^1   H_{t+1}^0    + N_{t+1}^0 H_{t+1}^1  \mid {\mathfrak A}_t \right)
%+ {\mathbb E} \left( N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right)  + L_t^1 \\
%& = {\mathbb E} \left[  N_{t+1}^0 \left(H_{t+1}^1  +  Q_{t+1}^1 H_{t+1}^0 \right) \mid {\mathfrak A}_t \right] + L_t^1 \\
%& = 0
%\end{align*}
%```
%where we used the implication that $ H_{t+1}^0  + L_{t+1}^0 = 0 $. The contribution:
%```{math}
%{\mathbb E} \left(N_{t+1}^0  H_{t+1}^1   \mid {\mathfrak A}_t \right) + L_t^1
%```
%is of the form used  for the first-order approximation without the recursive utility modification, except that the expectation is evaluated under the probability measure implied by $ N_{t+1}^0 $.  The recursive utility adjustment has us include the additional term:
%```{math}
%{\mathbb E} \left( N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right) 
%= {\frac {(\rho-1)}{2(1-\gamma_o)}} |\mu^o|^2H_{t+1}^0 \eqdef {\overline H}^1
%```
%which is constant over time.  Thus we write the first-order subsystem of equations as:
%```{math}
%{\mathbb E} \left(N_{t+1}^0   H_{t+1}^1 \mid {\mathfrak A}_t \right) + L_{t}^1 + {\overline H}^1 
% = 0.
%```
%We add to this the first-order approximation of the state dynamics inclusive of jump variables and evaluate expectations under the $ N_{t+1}^0 $ change of probability measure.  Thus the one-period conditional expectation of $ W_{t+1} $ is $ \mu^0 $. 
%
%(sec4)= 
%### Order two 
%
%The order two approximation of the product: $ N_{t+1} Q_{t+1} H_{t+1} + N_{t+1} L_{t+1} $ is:
%```{math}
%\begin{align*}
%& N_{t+1}^0 H_{t+1}^2 + N_{t+1}^2H_{t+1}^0  + 
%2 N_{t+1}^1 H_{t+1}^1 + L_t^2 \\
%+ 
%& 2 N_{t+1}^1  Q_{t+1}^1 H_{t+1}^0  +   2 N_{t+1}^0 Q_{t+1}^1 H_{t+1}^1 + N_{t+1}^0 Q_{t+1}^2 H_{t+1}^0 
%\end{align*}
%```
%The term $ N_{t+1}^2H_{t+1}^0 $ is zero and the term 
%```{math}
%{\mathbb E} \left( N_{t+1}^0 H_{t+1}^2  \mid {\mathfrak A}_t \right) + L_t^2
%```
%coincides with the contribution for the  second-order approximation  abstracting from recursive utility but evaluated under the change of measure induced by $ N_{t+1}^0 $.   
%Express $ H_{t+1}^1 $  as 
%```{math}
%:label: H1
%\begin{align} 
%H_{t+1}^1 & = \Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right). 
%\end{align}
%```
%We now consider the additional terms
%```{math}
%\begin{align} 
%2{\mathbb E}\left( N_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 
%(1 - \gamma_o)  {\mathbb E} \left[ N_{t+1}^0   \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right)H_{t+1}^1 \mid{\mathfrak A}_t \right], \\ 
%& = (1 - \gamma_o)  \Theta_2^1 \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2 \right) \\
%2 {\mathbb E} \left(N_{t+1}^1  Q_{t+1}^1 H_{t+1}^0 \mid {\mathfrak A}_t \right) & = (\rho  - 1)(1 - \gamma_o){\mathbb E}\left[ N_{t+1}^0 \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right) \left( {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1\right)  \mid {\mathfrak A}_t \right] H_{t+1}^0 \\
%& = (\rho  - 1) \mu^o \cdot \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2 \right) H_{t+1}^0 \\
%2 {\mathbb E} \left(N_{t+1}^0 Q_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 2(\rho - 1) {\mathbb E}\left[ N_{t+1}^0\left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right) \right] \mid {\mathfrak A}_t \right] \\
%&  = 2 \frac {(\rho - 1)}{(1-\gamma_o)}\left[   \Theta_2^1 \mu^0 + {\frac 1 2} \mu^o \cdot \mu^o \left(\Theta_0^1 + \Theta_1^1 X_t^1\right) \right] \\
%{\mathbb E} \left(  N_{t+1}^0 Q_{t+1}^2 H_{t+1}^0 \mid {\mathfrak A}_t \right) & =  (\rho - 1)^2{\mathbb E} \left[N_{t+1}^0 \left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2  \mid {\mathfrak A}_t \right] H_{t+1}^0 \\
%& = \left(\frac {1-\rho}{1-\gamma_o}\right)^2\left(|\mu^0|^2 + {\frac 1 4} |\mu^0|^4 \right) H_{t+1}^0
%\end{align} 
%```
%Denote the sum of the four terms in the second affine as ${\overline H}_t^2$. This random variable will be affine in 
%$ X_t^1 $, with a dynamic evolution determined by solving the first-order approximation.  Thus we write the subsystem of equations to be solved as:
%```{math}
%{\mathbb E} \left( N_{t+1}^0 H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\overline H}_t^2  = 0.
%```
%We add to this second-order subsystem, the second-order approximation of the state dynamics inclusive of the jump variables.  We substitute in the solution for the first-order approximation for the jump variables into both the first and second-order approximate state dynamics.  In solving the second-order jump variable adjustment we use expectations induced by $ N_{t+1}^0 $ zero throughout under which $ W_{t+1} $ is conditionally normally distributed with mean $ \mu^0 $ and covariance $ I $.
%
%
%### Steps for implementation 
%
%
%
%We implement these methods for second-order approximation using the following steps.
%
%
%
%1. Solve $H_{t+1}^0$ and $L_{t+1}^0$ for order zero state and jump variables. The outcome will be state invariant.
%
%2. Take as given a $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in representations {eq}`VminusR2`.
%
%3. Compute the first-order contribution to approximation by following the previous literature with expectations computed using the probabilities induced by $N_{t+1}^0$, which imply that $W_{t+1}$ has mean $\mu^0$. Express the solution as in {eq}`H1`.
%
%4. Compute the second-order contribution to the approximation by following the previous literature, again with the expectations induced by $N_{t+1}^0$.
%
%5. Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in representations {eq}`VminusR2` and return to {stepi}. Repeat until convergence.
%
%
%(sec2)= 
%## Appendix: Approximation formulas for expectation equations (approach two)
%
%For this solution, we iterate over $N_{t+1}^*$ approximation. Call the approximation ${\widetilde N}_{t+1}$ with an induced distribution for $W_{t+1}$ that is normal with conditional mean ${\tilde \mu}_t$ and covariance matrix ${\widetilde \Sigma}$. This distribution is used in both the first-order and second-order contributions to the approximation. The conditional mean for ${\tilde \mu}_t$ is affine in $X_t^1$. The following delineates the changes that need to be made.
%
%### First-order adjustment
%Compute:
%```{math}
%\begin{align*}
%{\mathbb E} \left( {\widetilde N}_{t+1}  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right) = \hspace{.2cm} &
%(\rho - 1) 
%{\mathbb E} \left[ {\widetilde N}_{t+1} \left({\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
%= \hspace{.2cm}&\left({\frac {\rho -1} {1-\gamma_o}}\right)  \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] H_{t+1}^0 \cr
% \eqdef \hspace{2cm} & {\widetilde  H}_t^1.
%\end{align*}
%```
%Then the equation to be solved is:
%```{math}
%{\mathbb E} \left({\widetilde N}_{t+1}   H_{t+1}^1 \mid {\mathfrak A}_t \right) + L_{t}^1 + {\widetilde H}^1_t  = 0.
%```
%
%### Second-order adjustment
%```{math}
%:label: second_affine_again
%\begin{align*}
% 2 {\mathbb E} \left({\widetilde N}_{t+1} Q_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) =  \hspace{.2cm} & 2(\rho - 1) {\mathbb E}\left[ {\widetilde N}_{t+1}\left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right) \right] \mid {\mathfrak A}_t \right] \cr
% = \hspace{.2cm}  & 2 \frac {(\rho - 1)}{(1-\gamma_o)}  \Theta_2^1 {\widetilde \Sigma} \mu^0 
%\cr  + &2 \frac {(\rho - 1)}{(1-\gamma_o)} \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1\left( {\tilde \mu}_t - \mu^0\right)  \right] \cr
%{\mathbb E} \left(  {\widetilde N}_{t+1} Q_{t+1}^2 H_{t+1}^0 \mid {\mathfrak A}_t \right) =  \hspace{.2cm}& (\rho - 1)^2{\mathbb E} \left[{\widetilde N}_{t+1} \left( 
%{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2  \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr 
%&  +   (\rho - 1) {\mathbb E} \left[{\widetilde N}_{t+1}  \left( 
%{\widehat V}_{t+1}^2  - {\widehat R}^2 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
%= \hspace{.2cm}& \left(\frac {1-\rho}{1-\gamma_o}\right)^2\left[{\mu^0}'{\widetilde \Sigma} \mu^0 + \left(\mu^0 \cdot {\tilde \mu}_t - {\frac 1 2} |\mu^0|^2 \right)^2\right]
%H_{t+1}^0 \cr 
% &  + \frac {(\rho - 1)} 2 \left[ \rm{tr}\left( \Upsilon_2^2{\widetilde \Sigma} - \Upsilon_2^2 \right)  
%  +  \left({\tilde \mu}_t - \mu^0\right)' \Upsilon_{2}^2  \left({\tilde \mu}_t - \mu^0\right)\right] H_{t+1}^0 \cr
%&  +  (\rho - 1) \left({\tilde \mu}_t  - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right) H_{t+1}^0.
%\end{align*}
%```
%
%Denote the sum of the two terms in {eq}`second_affine_again` as ${\widetilde H}_t^2.$ Then the equation to be solved is 
%```{math}
%{\mathbb E} \left( {\widetilde N}_{t+1} H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\widetilde H}_t^2  = 0.
%```
%
%### Updated recursive utility adjustments
%
%Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2,  \Upsilon_2^2$ used in representations {eq}`VminusR2`. Compute a new version of
%```{math}
%{\widetilde  N}_{t+1}  =  
%\frac 
%{\exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right]} 
%{{\mathbb E}\left( \exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right] \mid {\mathfrak A}_t \right)},
%```
%and deduce the implied ${\widetilde \mu}_t$ and ${\widetilde \Sigma}$. The conditional mean ${\tilde \mu}_t$ satisfies:
%```{math}
%{\widetilde \Sigma}^{-1} {\tilde \mu}_t =  \mu^0 +  {\frac {(1-\gamma_o)} 2} \left(\Upsilon_0^2 + \Upsilon_1^2 X_t^1 - \Upsilon_2^2 \mu^0 \right)
%```
%where the formula for ${\widetilde \Sigma}$ is
%```{math}
%\widetilde {\Sigma} =  \left[{\mathbb I} - {\frac {(1 - \gamma_o)}  2} \Upsilon_2^2\right]^{-1}.
%```
%
%With these adjustments, we iterate to convergence.
%

## Appendix D: Parameter values




To facilitate a comparison to a global solution method, we write down a discrete-time approximation to a continuous time version of such an economy. (See Section 4.4 of {cite}`KhorramiTourreHansen:2024` for a continuous-time benchmark model that our discrete-time system approximates.)  The parameter settings are:

|  $\eta_k$ | $\phi$  | $\beta_k$  | $\beta_1$ | $\beta_2 $ | $\mu_2$ |
| :------: | :-----: | :-----: | :-----: | :-----: | :-----: |
| .04 | 8  | .04 |  .056 |  .194 | $6.3\times10^{-6}$ |

\begin{align*}
\begin{bmatrix} \sigma_k \cr
\sigma_1 \cr
\sigma_2 
\end{bmatrix} = 
\sqrt{12} \begin{bmatrix}    .92 &   .40  &  0 \cr
 0 &  5.7 &  0  \cr
0 &  0 &  .00031 
\end{bmatrix} 
\end{align*}

The numbers for  $\eta_k, \phi, \beta_1,$ $\sigma_k$ and $\sigma_1$ are such that, when multiplied by stochastic volatility, they match the parameters from {cite}`LPH_TJS_tenuous`. In particular, the constant $Z^2$ which scales our $\sigma_k$ to match  is $7.6\times10^{-6},$ which is the 67th percentile of our $Z^2$ distribution. . While {cite}`LPH_TJS_tenuous` use a lower triangular representation for the two-by-two right block of 

```{math}
        \begin{bmatrix}
       \sigma_k \cr \sigma_1 \end{bmatrix} ,
```
we use an observationally equivalent upper triangular representation for most of the results. Finally, the numbers for $\beta_2$ and $\sigma_2$ come from {cite}`SchorfheideSongYaron:2018`, but they are adjusted for approximation purposes as described in Appendix A {cite}`KhorramiTourreHansen:2024`.    In both cases, we use the medians of their econometric evidence as input into our analysis.



























<!-- ```{bibliography}
:filter: docname in docnames
``` -->

