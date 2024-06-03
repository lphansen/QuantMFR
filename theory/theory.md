# Recursive Utility Framework

We aim to solve models which employ the following **recursive utility** framework:

$$
V_t = \left[ (1 - \beta) \left(C_t\right)^{1-\rho}  + \beta  \left( R_t \right)^{1-\rho} \right]^{\frac 1 {1-\rho}}
$$
where

$$\label{homog1b}
R_t = \left( {\mathbb E} \left[ \left( V_{t+1} \right)^{1-\gamma} \mid {\mathfrak A}_t \right] \right)^{\frac 1 {1-\gamma}} .
$$ 

Note that $0<\beta<1$ is a subjective discount factor, $\rho$ is the inverse elasticity of intertemporal substitution and $\gamma$ describes risk aversion. Since continuation values are determined only up to an increasing transformation, we work with logs of the above equation:

$$
\label{eq:value_recur5}
{\widehat  V}_t = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t] + \beta \exp \left[(1-\rho) {\widehat  R}_t \right] \right]
$$
where

$$
\label{eq:value_risk6}
{\widehat  R}_t = {\frac 1 {1 - \gamma}} \log {\mathbb E} \left(  \exp \left[ (1 - \gamma) {\widehat V}_{t+1} \right] \mid {\mathfrak A}_t \right).
$$

While $\gamma$ is usually interpreted as a measure of risk aversion, we now turn to an alternative interpretation which allows us to incorporate ambiguity aversion into our models.

## 1 Robustness

## 1.1 Motivation for Ambiguity Aversion

In early 2020, policymakers across the world were confronted with the decision of how to respond to the COVID-19 outbreak in the face of limited information regarding the disease's infectiousness and severity. Some countries opted for severe lockdowns, while others adopted a wait-and-see approach.

Similarly, designing the optimal climate policy requires reckoning with uncertainty surrounding the impact of economic activity on climate change as well as the adoption of new technologies. Should policymakers act under the assumption that a set of models are equally likely, or act under the worst-case scenario?

The answer depends on society's aversion to ambiguity and model misspecification.

## 1.2 Reinterpretation of Recursive Utility

.. (tbc)

Notice that $N_{t+1}^*$ has expectation one conditioned on ${\mathfrak A}_t$ and induces a change of probability measure motived by a robustness adjustment to potential model misspecification.    The expansion treats the approximation of $N_{t+1}^*$ and $ Q_{t+1}^*$ differently.  Both are computed in a separate step from the solution to equation (1) given $N_{t+1}^*$ and $Q_{t+1}^*$.  The implied one-period stochastic discount factor 

$$
\beta N_{t+1}^*Q_{t+1}^* \left( \frac {C_{t+1}}{C_t} \right)^{-\rho}
$$

where the last contribution is familiar from models with time separable, power utility preferences.  

We restrict $\gamma > 1$, and index it by the parameter ${\sf q}$ according to:

$$
\xi = \frac {1 - \gamma}{\sf q},  \;\;\;\; \xi_o = \frac {\sf q} {1 - \gamma_o} .  
$$

This embedding makes the uncertainty adjustments matter at a lower order.  It has a particularly nice interpretation when we interpret the recursion as a way to represent a preference for robustness to potential model misspecification.  We take ${\sf q} = 1$ as the economy of interest.


[1] Formally, ${ \frac {1}{\rho}}$ is the elasticity of intertemporal substitution.
<br><br>