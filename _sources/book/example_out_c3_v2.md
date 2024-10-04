(sec:chap2statincr)=
# Stationary Increments

Logarithms of many economic time series that appear to display stochastic growth can be modeled as having stationary increments. Multivariate versions of these models possess stochastic process versions of balanced growth paths. Applied econometricians seek permanent shocks that contribute to such growth. Furthermore, we shall see that it is convenient to pose central limit theory in terms of processes with stationary increments. The mathematical formulation in this chapter opens the door to studying these topics.
## Basic setup

We adopt assumptions from [](sec:inventing_past) that allow an infinite past and again let
$ {\mathfrak A}$ be a subsigma algebra of ${\mathfrak F}$ and 

```{math}
{\mathfrak A}_t = \left\{ \Lambda_t \in {\mathfrak F} : \Lambda_t = \{ \omega \in \Omega : {\mathbb S}^t(\omega)
\in \Lambda \} \textrm{ for some } \Lambda \in {\mathfrak F} \right\} .
```

Let $X$ be a scalar measurement function.
Assume that $Y_0$ is ${\mathfrak A}_0$ measurable and consider a scalar process $\{Y_t : t=0,1,... \}$ with stationary increments $\{X_t\}$:

```{math}
:label: eqn:add_1000
Y_{t+1} - Y_t = X_{t+1}
```

for $t=0,1, \ldots$. 
Let

```{math}
\nu = E\left(X_{t+1} \vert {\mathfrak I} \right),   \\
```
and 

```{math}
U_{t+1} = X_{t+1} - E\left(X_{t+1} \vert {\mathfrak A}_t \right).
```

We can interpret the above equations  as providing two contributions to the   $\{Y_{t}: t \ge 0\}$ process.
Thus, component $U_{t+1}$ is unpredictable and represents new information about $Y_{t+1}$ that arrives at date $t+1$.
Component $\nu$ is the trend rate of growth or decay in $\{Y_{t} : t \ge 0\}$ conditioned on the invariant events. In the following sections, we present a full  decomposition of a stationary increment process that will be useful both in connecting to sources of permanent versus transitory shocks and to central limit theorems.
## A martingale decomposition

A special class of stationary increment processes called additive martingales interests us.
````{prf:definition}
The process $\{Y_t^m : t=0,1,... \}$ is said to be an **additive martingale** relative to $\{ {\mathfrak A}_{t} : t=1,2,... \}$ if for $t=0,1,... $

- $Y_t^m$ is  ${\mathfrak A}_{t}$ measurable, and 
- $E\left(Y_{t+1}^m \vert {\mathfrak A}_t \right) = Y_t^m$ .
````

Notice that by the Law of Iterated Expectations, for a martingale  $\{Y_{t}^m :  t \ge 0\}$, best forecasts satisfy:
```{math}
E \left (Y_{t+j}^m \mid {\mathfrak A}_t \right) = Y_t^m
```
for $j \ge 1$. Under suitable additional restrictions on the increment process $\{X_t : t \ge 0 \}$, we can deploy a construction of {cite:t}`gordin` to show that the $\{V_t\}$ process contributes a martingale component to the $\{Y_t^m : t=0,1, ... \}$ process.[^gordinHallHeyde] Let ${\mathcal H}$ denote the set of all scalar random variables $X$ such that $E(X^2) < \infty$ and such that[^markovChapter]
```{math}
H_t = 
\sum_{j=0}^\infty E\bigl( X_{t+j} - \nu \vert {\mathfrak A}_t \bigr)
```
is well defined as a mean-square convergent series. Convergence of the infinite sum on the right side limits temporal dependence of the process $\{ X_t \}$. For example, it can exclude so-called long memory processes.[^longMemory]

Construct the one-period ahead forecast of  $H_{t+1}$: 
```{math}
H_t^+ = E\left( H_{t+1} \mid {\mathfrak A}_{t} \right)
``` 
Notice that
```{math}
X_t - \nu =  H_t - H_t^+ =  G_t  + \left( H_{t-1}^+ - H_t^+ \right)
```
where 
```{math}
:label: martdiff
G_{t} = H_{t} - H_{t-1}^+ = H_t - E\left( H_{t} \mid {\mathfrak A}_{t-1} \right). 
```
Since  $G_t$ is a forecast error,
```{math}
E \left( G_{t+1} \vert {\mathfrak A}_{t} \right) = 0.
```

Assembling these parts, we have
```{math}
:label: eqn:add1001
Y_{t+1} - Y_t = X_{t+1} =   \nu + G_{t+1}  + H_t^+ - H_{t+1}^+ .
```
Let
```{math}
Y^m_t = \sum_{j=1}^t G_j  .
```
Since $Y_t^m$ is ${\mathfrak A}_{t}$ measurable, the equality
```{math}
E \left(  \sum_{j=1}^{t+1}   G_j    \mid {\mathfrak A}_t \right) = \sum_{j=1}^{t}  G_j 
```
implies that the process $\{Y_t^m : t \ge 0 \}$
is an *additive martingale*.  

For a given stationary increment process, $\{Y_t : t \ge 0\}$, express the martingale increment as 
```{math}
:label: limitnews
G_{t}  = \sum_{j=0}^\infty \left[ E\left( X_{t+j}   \mid  {\mathfrak A}_{t} \right) -  E\left( X_{t+j}   \mid  {\mathfrak A}_{t-1}   \right) \right] 
= \lim_{j \rightarrow \infty} \left[ E\left(Y_{t+j} \vert {\mathfrak A}_{t} \right)  -  E\left(Y_{t+j} \vert {\mathfrak A}_{t-1} \right) \right] .
```
So the increment to the martingale component of $\{Y_t : t \ge 0 \}$ is new information about the limiting optimal forecast of $Y_{t+j}$
as $j \rightarrow + \infty$.

By accumulating equation {eq}`eqn:add1001` forward, we arrive at:

[^gordinHallHeyde]: Also see {cite:t}`hallheyde`.
[^markovChapter]: The random variable $H_t$ somewhat resembles an "undiscounted" version of the resolvent operator that plays an important role in the analysis of Markov processes in chapter [](chap:markov).
[^longMemory]: See, for instance, {cite:t}`GrangerJoyeux:1980`,  {cite:t}`GewekePorter-Hudak:1983` and {cite:t}`Robinson:1994`.

````{prf:proposition}
:label: decomp00
If $X$ is in ${\mathcal H}$, the stationary increments process $\{Y_t : t=0,1,...\}$ satisfies the additive decomposition

```{math}
\begin{matrix}
Y_{t} & = & \underbrace{t\nu} & + & Y_t^m & - &\underbrace{ H_t^+} & + & \underbrace{Y_0 + H_0^+}.\cr
&&\textbf{trend} &&\textbf{martingale}&& \textbf{stationary} && \textbf{invariant}
\end{matrix}
```
The stationary increment process, $\{Y_{t}^m : t\ge 0 \},$ is the martingale component with  $Y_0^m = 0$,   The component $\{H_{t}^+\}$ is stationary.  The other components are constant over time.
````

{prf:ref}`decomp00` decomposes  a stationary-increment process into a linear  time trend, a martingale, and a transitory component. A permanent shock is the increment to the martingale.  The martingale and transitory contributions are typically correlated.  
````{prf:example}
:label: ex:chap2ma

(Moving-average increment process) Consider again the {prf:ref}`ex:MA` moving-average process:

```{math}
:label: eqn:firstoderex
X_{t} = \sum_{j=0}^\infty \alpha_j \cdot W_{t-j} .
```

Use this $\{X_t\}$ process as the increment for $\{ Y_t : t \ge 0  \}$ in formula {eq}`eqn:add_1000`.  New information about the  unpredictable component of $X_{t+j}$ for $j \ge 0$  that arrives at date $t$ is 

```{math}
E \left( X_{t+j} \mid {\mathfrak A}_{t} \right)-  E \left( X_{t+j} \mid {\mathfrak A}_{t-1} \right)= \alpha_{j} \cdot W_{t}
```

Summing these terms over $j$ gives

```{math}
G_{t} = \alpha(1)  \cdot W_{t}  
```

where

```{math}
\alpha(1) = \sum_{j=0}^\infty \alpha_j
```

provided that the coefficient sequence $\{ \alpha_j : j\ge 0\}$
is  summable, a condition that restricts  temporal dependence of the increment process $\{X_t\}$.   Indeed, it is possible for $\alpha(1) = \infty$ or for it not to be well defined while

```{math}
 \sum_{j=0}^\infty |\alpha_j|^2 < \infty
```

ensuring that $X_t$ is well defined.  This possibility opened the door to the literature on long-memory processes that allow for $\alpha(1)$ to be infinite as discussed in {cite:t}`GrangerJoyeux:1980` and elsewhere.

In what follows, we presume that $\alpha(1)$ is finite.  This sum of the coefficients  $\{\alpha_j: j\ge 0 \}$ in  moving-average representation {eq}`eqn:firstoderex` for
the **first difference** $Y_{t+1}  - Y_t = X_{t+1}$ of $\{ Y_t : t=0,1,.... \}$ tells the **permanent** effect of $W_{t+1}$ on current and future values of the
**level** of $Y$, i.e., the effect on $\lim_{j\rightarrow + \infty} Y_{t+j}$. 
Models of {cite:t}`blanchardquah` and {cite:t}`shapirowatson` build on this
property.

The variance of the random variable $\alpha(1) \cdot W_{t+1}$  conditioned on the invariant events in ${\mathfrak I}$  is $|\alpha(1)|^2$.  The overall variance of $X_{t}$ is given by

```{math}
\sum_{j=0}^\infty|\alpha_j|^2 \ne |\alpha(1)|^2.
```

To form a permanent-transitory shock decomposition, construct the  permanent shock as:

```{math}
W_{t+1}^p = \left( \frac {1}{|\alpha(1)|} \right) \alpha(1) \cdot W_{t+1}
```
where we introduce an additional scaling so the permanent shock has variance one.  Form

```{math}
W_{t+1}^{tr} = W_{t+1} -   \left( \frac {1}{|\alpha(1)|} \right) \alpha(1) W_{t+1}^p 
```
which by construction will be uncorrelated with $W_{t+1}^p$.  Since  the covariance matrix of $W_{t+1}^{tr}$ will be singular, the components of $W_{t+1}^{tr}$ can be expressed as linear combinations of a vector of transitory shocks with unit variances.
````


````{prf:example}
:label: ex:helpBQ

This is a process in which $W_t$ has transient but no permanent effects on future $Y$'s. Let $\alpha_0 = 1$ and $\alpha_j = (\lambda -1) \lambda^{j-1}$ for $j \geq 1$ and $-1 < \lambda < 1$. Construct the power series

```{math}
:label: factor
\alpha(\zeta) = 1 - \sum_{j=1}^\infty  ( 1- \lambda ) \lambda^{j-1} \zeta^j = 1 - {\frac { (1 - \lambda) \zeta}{ 1 - \lambda \zeta}} = {\frac {1 - \zeta}{1 - \lambda \zeta}} .
```

Evidently, $\alpha(1) = 0$. Define

```{math}
H_t  = - \sum_{j=0}^\infty \lambda^j W_{t-j}
```

and note that since {eq}`factor` is satisfied

```{math}
Y_{t+1} - Y_t = -H_{t+1} + H_t.
```

The process $\{ Y_t : t = 0,1,... \}$ is stationary provided that $Y_0 = - H_0$, ensuring that $Y_t = - H_t$ for all $t \ge 0$.
````
## Central limit approximation

{prf:ref}`ex:chap2ma` starts from a moving average of martingale differences that is used as an increment $\{X_t \}$
to a $\{Y_t: t \ge 0\}$ process, after which it constructs a process of innovations
to the martingale component of the $\{Y_t: t \ge 0  \}$ process. That analysis illustrates the workings of
an
operator $\mathbb{D}$ that maps an admissible increment process in $\mathcal{H}$ into the innovation
in a martingale component. To construct $\mathbb{D}$, let $\mathcal{G}$ be the set of all random variables $G$ with finite second moments that satisfy the conditions that
$G$ is $\mathfrak{A}$ measurable and that $E(G_1 \vert \mathfrak{A}) = 0$ where $G_1 = G \circ \mathbb{S}$. Define $\mathbb{D}: \mathcal{H} \rightarrow \mathcal{G}$ via
```{math}
\mathbb{D}(X) = G .
```
Both $\mathcal{G}$ and $\mathcal{H}$ are linear spaces of random variables and $\mathbb{D}$ is a linear transformation.
The operator $\mathbb{D}$ plays a prominent role in a central limit approximation.

To form a central limit approximation, construct the following scaled partial sum that nets out trend growth
```{math}
{\frac 1 {\sqrt{t}}}(Y_t - \nu t) = {\frac 1 {\sqrt{t}}} Y_t^m - {\frac 1 {\sqrt t}} H_{t}^+
+ {\frac 1 {\sqrt{t}}} (H_0^+ + Y_0)  
```
where
```{math}
Y_t^m= \sum_{j=1}^t G_j
```

From {cite:t}`billingsley`'s central limit theorem for martingales
```{math}
{\frac 1 {\sqrt{t}}} Y_t^m \Rightarrow \mathcal{N} \left(0, E\left[ \mathbb{D}(X)^2 \vert \mathfrak{I} \right] \right)
```
where $\Rightarrow$ denotes weak convergence, meaning convergence in distribution. Clearly, $\{(1/ {\sqrt t}) H_{t}^+\}$ and $\{(1/{\sqrt{t}}) (H_0^+ + Y_0) \}$ both converge in mean square to zero.
````{prf:proposition}
:label: prop:gordin
For all stationary increment processes $Y_t : t=0,1,2, ...$ represented by $X$ in $\mathcal{H}$

```{math}
{\frac 1 {\sqrt{t}}}(Y_t - \nu t) \Rightarrow  {\mathcal{N}} \left( 0,
E\left[ \mathbb{D}(X)^2 \vert {\mathfrak{I}} \right] \right) .
```

Furthermore,

```{math}
E\left[ \mathbb{D}(X)^2 \vert {\mathfrak{I}} \right]
= \lim_{t \rightarrow \infty}  E \left[ \left({\frac 1 {\sqrt{t}}}  \left(Y_t - t \nu \right)  \right)^2 \Bigl| {\mathfrak{I}} \right].
```
````
## Cointegration

Linear combinations of stationary increment processes $Y_t^1$ and $Y_t^2$ have stationary increments. For real-valued scalars $r_1$ and $r_2$, form
```{math}
Y_{t} = r_1 Y_{t}^1 + r_2 Y_{t}^2
```
where
```{math}
\begin{align*}
Y_{t+1}^1 - Y_t^1 & = X_{t+1}^1 \\
Y_{t+1}^2 - Y_t^2 & = X_{t+1}^2.
\end{align*}
```
The increment in $\{Y_t : t=0, 1, \ldots \}$ is
```{math}
X_{t+1} = r_1 X_{t+1}^1 + r_2 X_{t+1}^2
```
and
```{math}
Y_0 =  r_1 Y_0^1 + r_2 Y_0^2.
```
The {prf:ref}`decomp00` martingale component of $\{ Y_t : t \ge 0 \}$ is the corresponding linear combination of the martingale
components of $\{ Y_t^1 : t =0,1,...\}$ and $\{ Y_t^2 : t =0,1,...\}$. The {prf:ref}`decomp00`  trend component of $\{ Y_t : t =0,1,  \ldots \}$ is the corresponding linear combination of the
trend  components of $\{ Y_t^1 : t =0,1, \ldots \}$ and $\{ Y_t^2 : t =0,1,  \ldots \}$.

{prf:ref}`decomp00` sheds light on the cointegration concept of {cite:t}`englegranger` that is associated with linear combinations of stationary increment processes whose trend and martingale components are both zero. {citeauthor}`englegranger` Call two processes *cointegrated* if there exists a linear combination of them that is stationary.[^cointegration-def] That situation prevails when there exist real-valued scalars $r_1$ and $r_2$ such that
```{math}
\begin{eqnarray*}
r_1 \nu_1 + r_2 \nu_2 & = & 0 \\
r_1 \mathbb{D}(X^1) + r_2 \mathbb{D}(X^2) & = & 0,
\end{eqnarray*}
```
where the $\nu$'s correspond to the trend components in {prf:ref}`decomp00`. These two zero restrictions imply that the time trend and the martingale component for the linear combination $Y_t$ are both zero.[^cointegration-vector] When $r_1 = 1$ and $r_2 = - 1$, the stationary increment processes $Y_{t}^1$ and $Y_{t}^2$ share a common growth component.

This notion of cointegration provides one way to formalize balanced growth paths in stochastic environments through determining a linear combination of growing time series for which stochastic growth is absent.

[^cointegration-def]: The {cite:t}`BoxTiao:1977` ‘‘canonical correlation’’ approach to linear time series analysis anticipated, at least partially, the co-integration restrictions of time series econometricians and macroeconomists.
[^cointegration-vector]: The  cointegration vector $(r_1, r_2)$ is  determined  only up to scale.

<!-- 
## Permanent Income Models

This section describes a class of consumption-savings models that we cast in terms of exogenous and endogenous variables all of which are additive processes with stationary increments like {eq}`eqn:add_1000`. These models generalize one with which we worked in {cite:t}`HSRoberds`, {cite:t}`robustpi`, and {cite:t}`HansenSargent_Recursive_Models[ch.~11]`. Those studies had formulated versions of the models to be presented here but with exogenous and endogenous processes that were constrained to be stationary Markov processes.

In the models here,
restrictions across exogenous and endogenous variables emerge from two theoretical sources: (1) present value budget balance, and (2) an optimizing consumer's first-order necessary condition that implies that marginal utilities form a martingale.

### Present-value budget balance

Let $Y_t$ be time $t$ nonfinancial income minus time $t$ expenditures and assume that it is a process with stationary increments:

```{math}
:label: eqn:permincome101
Y_{t+1} - Y_t = X_{t+1} + \eta,
```

where

```{math}
:label: eqn:permincome102
X_{t+1} = \sum_{j=0}^\infty \alpha_j \cdot W_{t+1-j}.
```

We call $Y_t$ the flow *surplus* and $- Y_t$ the flow *deficit*. Suppose that $\lambda \in (0,1)$ is a constant discount factor and assume that a risk-free asset has a gross return $\lambda^{-1}$.[^footnotelambda]

Suppose that the riskless asset is a consumer's only vehicle for borrowing or lending, so that

```{math}
:label: eqn:pvdiff
\lambda \left(Y_{t+1} + K_{t+1} \right) = K_t ,
```

where $K_t$ is a stock of the asset at the end of period $t$.

It can be verified that equation {eq}`eqn:pvdiff` implies

$$
\lambda \left( Y_{t+1} - Y_t + K_{t+1} - \chi Y_{t} \right) = K_t - \chi Y_t ,
$$

or equivalently

$$
\chi \left( Y_{t+1} - Y_t \right) + \lambda \left(K_{t+1} - \chi Y_{t+1} \right) = K_t - \chi Y_t.
$$

where $ \chi = \lambda + \lambda \chi $.

Since $\lambda \in (0,1)$, $ \chi > 0$. We can solve this first-order linear difference equation in $\{K_t - \chi Y_t\}$ forward to get

```{math}
:label: eqn:PVBB201
K_t - \chi Y_t  = \chi \sum_{j=0}^\infty \lambda^j (Y_{t+j+1} - Y_{t+j} ) ,
```

an equation that we require to hold for every sample path of $\{Y_t\}$. We regard equation {eq}`eqn:pvdiff` as a time $t$ flow budget constraint and equation {eq}`eqn:PVBB201` as a time $t$ expression of a present value budget constraint. Because equation {eq}`eqn:PVBB201` holds for every sample path, it is legitimate to average over all such paths with an arbitrary probability distribution. By averaging with the probability distribution governing the $\{Y_t\}$ process conditional on ${\mathfrak A}_{t+1}$, we deduce

```{math}
:label: eqn:PVBB202
K_t - \chi Y_t = \chi \left[ E \sum_{j=0}^\infty \lambda^j \left(Y_{t+j+1}  - Y_{t+j}\right) \Bigr| {\mathfrak A}_{t+1} \right].
```

Equation {eq}`eqn:PVBB202` imposes restrictions that we summarize as follows.[^footnoteHS1980]

[^footnotelambda]: In section XXXXXX, we will relax constancy of $\lambda$.

[^footnoteHS1980]: {cite:t}`HS1980` employ a closely related frequency-domain argument.
````{prf:theorem
:label: prop:alpha0
Representation {eq}`eqn:PVBB201` implies the moving average representation
```{math}
K_t - \chi Y_t = \sum_{j=0}^\infty  \alpha^*_j W_{t+1-j} + \frac{\chi}{1-\lambda} \eta
```
where the $z$-transform of the $\{\alpha^*_j\}_{j=0}^\infty$ sequence satisfies
```{math}
\alpha^*(\zeta) =  \frac {\chi \zeta }{(1-\lambda) (\zeta - \lambda)}\alpha(\zeta) -  \frac  { \chi\lambda }{(1-\lambda)(\zeta - \lambda)}\alpha(\lambda),
```
where $\zeta$ is a complex number. The complex-valued function $\alpha^*(\zeta)$ of $\zeta$ has a removable singularity at $\zeta = \lambda$. Since $K_t - \chi Y_t$ depends on date $t$ information only and not on the shock at date $t+1$, $\alpha^*(0) = 0$. It follows that $\alpha(\lambda)= 0$ and that consequently
```{math}
\alpha^*(\zeta) = \frac {\chi \zeta }{(1-\lambda)(\zeta - \lambda)}\alpha(\zeta) .
```
````

````{prf:theorem}
:type: remark

The Proposition {prf:ref}`prop:alpha0` restriction that $\alpha(\lambda) = 0$ is an implication of intertemporal budget constraint {eq}`eqn:PVBB202` that restricts the expected present value of the surplus process $\{Y_t\}$. The following argument cast in the time domain substantiates this interpretation. A shock $W_{t+1}$ at date $t+1$ contributes

```{math}
\frac{1}{1 - \lambda} \lambda^j \alpha_j \cdot W_{t+1}
```

to the deficit at date $t+j+1$. Discounting and summing across future dates should give an outcome that equals zero for almost all shock realizations:

```{math}
:label: eqn:PVBB101
0 = \frac{1}{1 - \lambda} \sum_{j=0}^\infty \lambda^j \alpha_j = \frac{1}{1-\lambda} \alpha(\lambda).
```

*Note: Lars XXXXX: do we want to transform the "for almost all" statement into a statement about "almost everywhere" with respect to some probability measure?* This is another way of confirming that $\alpha(\lambda) = 0$.
````
````{prf:remark}
The right side of equation {eq}`eqn:PVBB201` is stationary, implying that $(1 , - \chi)$ is a cointegrating vector for the $\{K_t, Y_t\}$ process.
````

### Additional model components

A present value budget balance restriction like equation {eq}`eqn:PVBB101` is an important component of a class of permanent income model. Other model components determine the consumption expenditures that contribute negatively to the surplus. To construct such a model, we take as exogenous inputs

```{math}
:label: eqn:statinc1
L_{t+1} - L_t  = X_{t+1}^{[L]} + \eta_1\\
B_{t+1} - B_t  = X_{t+1}^{[B]} + \eta_2 ,
```

where $\eta_1, \eta_2$ are scalar constants, $L_t$ is income at $t$, $B_t$ is a preference shock, and

```{math}
:label: eqn:statinc2
X_{t+1}^{[L]}  = \sum_{j=0}^\infty  \alpha^{[L]}_j   \cdot W_{t+1-j} \\
X_{t+1}^{[B]}  = \sum_{j=0}^\infty \alpha^{[B]}_j  \cdot W_{t+1-j} .
```

To allow for habit persistence and durable goods, we assume that a "consumption services" process $\{S_t\}$ is governed by

```{math}
S_t  = (1 + \kappa) C_t - \kappa H_{t-1}
```
{:label: service}

```{math}
H_t  = \exp(-\delta_h) H_{t-1} + [ 1 - \exp(-\delta_h)]  C_t
```
{:label: householdstock}

where $\delta_h  \geq 0$ is a depreciation rate, so that $H_t$ is a geometrically weighted average of current and past consumption $C_t$. Setting $\kappa > 0$ makes consumption services depend positively on current consumption and negatively on past consumption, a pattern of dependence symptomatic of "habit persistence." Setting $\kappa \in (-1,0)$ makes consumption services depend positively on both current and past consumption, indicating "durable" consumption goods.

We also assume that there is a linear asset accumulation technology

```{math}
C_t + K_{t} - K_{t-1} = \rho K_{t-1} + L_t,
```

where $K_t$ is asset holdings at the end of period $t$. We can rewrite the preceding equation as

```{math}
K_t = (1 + \rho) K_{t-1} + L_t - C_t .
```

Solving forward gives

```{math}
K_{t-1} = \sum_{j=0}^\infty \lambda^{j+1} \left(  C_{t+1} -L_{t+j} \right) ,
```

where $\lambda = \frac{1}{1 + \rho}$. This equation implies the following version of the present-value budget balance restriction {eq}`eqn:PVBB101` for the net surplus process  $Y_t = L_t - C_t$, $t=0,1,...$:

```{math}
\alpha^{[L]}(\lambda) - \alpha^{[C]}(\lambda)  =0 .
```

(sec:MUSmart)=
### Marginal utility of services is a martingale

At this point, we jump ahead and complete our model by just arbitrarily restricting a process that describes the marginal utility of services experienced by an optimizing consumer. In subsection [sec:MUprocessopt](sec:MUprocessopt) we shall describe a consumer's optimum problem that rationalizes this process.

Let the marginal utility of consumption services be:

```{math}
MS_t =   \left( B_t - S_t \right),
```

so that the process  $B_t$ shifts the marginal utility for consumption services. Temporarily suppose that the marginal utility process is a martingale[^footMartingale]:

```{math}
:label: eqn:MSmart100
MS_{t+1} - MS_t = \alpha^{[M]} \cdot W_{t+1}.
```

From equation {eq}`service` and  $ MS_t = \left( B_t - S_t \right)$, we deduce

```{math}
:label: eqn:MSrep100
MS_t   = B_t - (1 + \kappa) C_t + \kappa H_{t-1} .
```

We seek a consumption process that is an additive process driven by stationary increments and so has the form

```{math}
:label: eqn:Caddprocess1
C_{t+1} - C_t = \sum_{j=0}^\infty \alpha_j^{[C]} W_{t+1 -j} .
```

The first difference of equation {eq}`eqn:MSrep100` implies that the following equation connects $\zeta$-transforms of the increments to $MS, B, C,$ and $H$:

```{math}
(1 + \kappa) \alpha^{[C]}(\zeta) - \zeta\kappa \alpha^{[H]}(\zeta)  =  \alpha^{[B]}(\zeta) - \alpha^{[M]} 
```

Equation {eq}`householdstock` implies the following restriction on $\zeta$-transforms of increments to $H$ and $C$:

```{math}
\left[1 -  \exp(-\delta_h) \zeta \right] \alpha^{[H]}(\zeta)  - \left[1 - \exp(-\delta_h) \right] \alpha^{[C]}(\zeta) = 0.
```

Combining the two preceding equations, we assemble the following system:

```{math}
\begin{bmatrix}
 (1 + \kappa)  & - \zeta\kappa  \\
- [1 - \exp(-\delta_h)] &  [1 -  \exp(-\delta_h) \zeta ] \end{bmatrix}
\begin{bmatrix} \alpha^{[C]}(\zeta) \\ \alpha^{[H]}(\zeta) \end{bmatrix}
 = \begin{bmatrix} \alpha^{[B]}(\zeta) - \alpha^{[M]} \\ 0 \end{bmatrix} . 
```

(sec:MUprocessopt)=
### Martingale marginal utility process

We now pose an optimum problem that implies that $\{ MS_t : t=0,1,...\}$ is a martingale as we assumed in subsection [Marginal utility of services is a martingale](sec:MUSmart). Let the optimal value function $V_t$ for a representative consumer satisfy the Bellman equation

```{math}
:label: eqn:Bellman_MUmart
V_t = - {\frac {1 - \exp(-\delta)} 2} MS_t \cdot MS_t + \exp(-\delta) {\mathcal R}_t(V_{t+1})
```

where the risk sensitivity operator ${\mathcal R}_t$ is defined as

```{math}
{\mathcal R}_t(V_{t+1}) =  - {\frac 1 {1 -\gamma}} \log E\left[ \exp\left(  (1 - \gamma) V_{t+1}\right)  \vert {\mathfrak A}_t \right]
```

and

```{math}
MS_{t+1} = MS_t + \alpha^{[M]} \cdot W_{t+1},
```

and where $W_{t+1}$ is described by a multivariate normal distribution with mean zero and an identity as the covariance conditioned on ${\mathfrak A}_t$. The risk-sensitive continuation value function ${\mathcal R}_t(V_{t+1})$ is itself an indirect utility function from the following minimization problem:

```{math}
{\mathcal R}_t(V_{t+1}) =
\min_{N_{t+1}}  E\left(N_{t+1} V_{t+1}  \vert {\mathfrak A}_t \right) + {\frac 1 {1 - \gamma}} E\left( N_{t+1} \log N_{t+1} \vert {\mathfrak A}_t \right) ,
```

where the minimization is subject to $N_{t+1} \ge 0$, $E\left( N_{t+1} \vert {\mathfrak A}_t \right) = 1$. We can interpret $N_{t+1}$ as a likelihood ratio that alters the conditional distribution of the shocks $W_{t+1}$. The minimizing  $N_{t+1}$ is

```{math}
N_{t+1} = {\frac  {\exp\left[(1 - \gamma) V_{t+1} \right]} { E \left( \exp\left[(1 - \gamma) V_{t+1} \right] \vert {\mathfrak A}_t \right)}} ,
```

#### Euler equations

The marginal utility of consumption $C_t$ is proportional to the marginal utility of services, so it is also a martingale when the marginal utility of consumption services $S_t$ is a martingale. Note that

```{math}
E\left( N_{t+1} MS_{t+1} \vert {\mathfrak A}_t \right)  =
\left[ 1 + {\sf s} (\gamma - 1) \left(\alpha^{[M]}\right)'\left[I + (1-\gamma){\sf s} \alpha^{[M]}{\alpha^{[M]}}'\right]^{-1} \alpha^{[M]}\right] MS_t .
```

Since the marginal utility for consumption is proportional to $MS_t$ and since the gross marginal product of capital is $1 + \rho$, the Euler equation for accumulating capital is

```{math}
\exp(-\delta)(1 + \rho)  E\left( N_{t+1} MS_{t+1} \vert {\mathfrak A}_t \right) = MS_t .
```

### Stochastic discounting

{cite:t}`robustpi` impute a robust savings problem like the one constructed in subsection [Martingale marginal utility process](sec:MUprocessopt) to a planner in a representative agent economy. The planner's aversion to model misspecification affects equilibrium valuations of risky payoff streams in ways that we shall describe in chapter XXXXXX.

[^footMartingale]: {cite:t}`HansenSargent_Recursive_Models` present models in which such a marginal utility process is a martingale as we assume here.
## Long-term consumption risk

{cite:t}`hhl` had the idea of using covariation with other time series to help infer long-run stochastic components of consumption.
Figure plots logarithms of nondurable consumption $C_t$ and corporate earnings $N_t$. The absence of an obvious trend or martingale in the second panel, which plots the difference between the logarithms of nondurable consumption and corporate earnings, suggests the presence of common trend and martingale components in the two series themselves, an observation that led {cite:t}`hhl` to impose co-integration between the logarithms of consumption and corporate earnings and thereby restrict them to grow together. $\{W_{t+1}\}$ is an i.i.d. sequence of ${\cal N}(0,I)$ random vectors; then to choose $X_t$ to have the growth rate of consumption (expressed in logarithms) as its first entry and the logarithm of corporate earnings minus the logarithm of consumption in the second position, then to fill in the remaining components of $X_t$ with lags of these and any other variables that help forecast the logarithms of corporate earnings and consumption.
This specification leaves us with two additive functionals with increments:

```{math}
\log Y^{[1]}_{t+1} - \log Y^{[1]}_t  = \nu_1 + X_{t+1}^{[1]} 
```
```{math}
\log Y^{[2]}_{t+1} - \log Y^{[2]}_t  = \nu_2  + X_{t+1}^{[2]} - X_t^{[2]}  + X_{t+1}^{[1]}
```

where $Y^{[1]}_{t+1} = C_{t+1}$ and $ Y^{[2]}_{t+1} = N_{t+1}$.
 The two additive functionals $\{\log Y^{[1]}_{t+1}\}$ and $\{\log Y^{[2]}_{t+1}\}$  share the same martingale and trend components but have different transitory components.

Notice that

```{math}
\log Y^{[1]}_{t+1} - \log Y^{[1]}_t  = \nu_1 + D \cdot X_t + F \cdot W_{t+1}
```
where
```{math}
D = A'U_1,  \  \  F = B' U_1
```
and $U_1$ is a vector of zeros except for a one in the first position. The impulse response vector of $\log C_{t+1} - \log C_t$ to $W_{t+1}$ is

```{math}
F,  B'D,  B'A'D, \cdots ,
```
which has $z$-transform:
```{math}
F + \zeta B'(I - \zeta A')^{-1} D = B' (I - \zeta A')^{-1} U_1.
```
The martingale increment scaled to have unit standard deviation is $F^* \cdot W_{t+1}$, where
```{math}
F^* = {\frac 1 {|B'(I -  A')^{-1} U_1 |}}  B' (I -  A')^{-1} U_1
```
and the $z$-transform of the impulse response function of $\log C_{t+1}$ to the martingale increment is
```{math}
\left( {\frac 1 {1 - \zeta}}\right) F^* \cdot [B' (I - \zeta A')^{-1} U_1].
```
(sec:Camp_Shiller)=
## Loglinear Approximations

It is common to use log-linear approximations to describe how payouts and required return processes affect asset values.
Let $V_{t}$ be the value at date $t$ of a claim to a stream $\{G_{t+j}\}_{j=0}^\infty$ of positive cash flows.
A one-period gross return on the asset is

```{math}
R_{t+1} = {\frac {V_{t+1} + G_{t+1}} {V_t}} =  \left({\frac {V_{t+1}/G_{t+1}  + 1} {V_t/G_t}}\right) \left({\frac {G_{t+1}}{G_t}} \right).
```

Let

```{math}
\lambda = {\frac 1 {\exp(\mu_v) + 1}} .
```

Approximate $\log V_t - \log G_t$ around its mean  $\mu_v$ to obtain

```{math}
\log \left(V_{t+1}/G_{t+1} + 1 \right) \approx
\log \left[ \exp\left( \mu_v \right) + 1 \right] +
 \lambda \left( \log V_{t+1} - \log G_{t+1} - \mu_v \right)  .
```

Take the implied difference equation

```{math}
\log V_t  - \log G_t & = \log \left[ \exp\left( \mu_v \right) + 1 \right]  + \cr &
 \lambda \left( \log V_{t+1} - \log G_{t+1} - \mu_v \right)   + \log G_{t+1} - \log G_t - \log R_{t+1}
```

and solve it forward to obtain the "present value model"

```{math}
:label: eqn:CSpresentval
\log V_t - \log G_t - \mu_v  & = \log \left[ \exp\left( \mu_v \right) + 1 \right] - \mu_v \cr &   + \sum_{j=0}^\infty \lambda^j
\left( \log G_{t+1+j} - \log G_{t+j} - \log R_{t+1+j} \right) .
```

Let $\mu_r$ be the mean of  $\log R_{t+1}$ and let  $\mu_g$ be the mean of
$\log G_{t+1} - \log G_t$.  Then

```{math}
- \mu_v = \log \left[ \exp\left( \mu_v \right) + 1 \right] - \mu_v - {\frac {\mu_g - \mu_r} {1 - \lambda}}
```

or

```{math}
\log \left[ \exp\left( \mu_v \right) + 1 \right] = {\frac {\mu_g - \mu_r} {1 - \lambda}},
```

which expresses  $\mu_v$ implicitly as a nonlinear function of $\mu_g$ and $\mu_r$.

To turn these calculations into a complete model of a $\{V_t, G_t, R_t\}$ process, we must specify stochastic processes for $\{\log R_t\}$ and $\{G_t\}$, the  "forcing" processes  appearing on the right side of equation {eq}`eqn:CSpresentval`. Here we assume the following two moving-average representations:

```{math}
:label: eqn:asset_price_inputs
\log R_{t+1} & =  \mu_r + \alpha^r({\mathcal L}) W_{t+1} \cr
\log G_{t+1} - \log G_t & =  \mu_g + \alpha^g({\mathcal L}) W_{t+1}
```

where ${\mathcal L}$ is again the lag operator.
We seek a moving average representation for $\log V_t - \log G_t$ implied by the processes {eq}`eqn:asset_price_inputs`.
We guess that it has the form

```{math}
:label: eqn:asset_price_outcome
\log V_{t} - \log G_t  = \mu_v +  \alpha^v({\mathcal L}) W_t. 
```

To fill in objects in {eq}`eqn:asset_price_inputs`, we assume that $\log R_t$ is a stationary process and that $\log G_t$ is a process with stationary increments. In writing form {eq}`eqn:CSpresentval`, we are guessing that $\log V_t - \log G_t$ is a stationary process and that $\log V_t$ is a process with stationary increments that is cointegrated with $\log G_t`.{^footnotelabel}` See [Cointegration](sec:cointegrate) for a definition  of cointegration.
The present-value model {eq}`eqn:CSpresentval` obtained from the log-linear approximation implies the cross-equation restriction that we again express in terms of $\zeta$-transforms

```{math}
\alpha_v(\zeta) = \left({\frac 1 {\zeta - \lambda}}\right) \left[ \alpha_g(\zeta) - \alpha_r(\zeta)\right] 
```

where $\zeta$ is a complex-valued scalar.
The assumption that  $\log V_t - \log G_t$  depends only on current information implies that $\alpha_v$ has a power series representation valid for $| \zeta | < 1$, so that  $\zeta = \lambda$ is a removable singularity, which in turn implies that

```{math}
\alpha_g(\lambda) = \alpha_r(\lambda).
```

This equation asserts that, when discounted appropriately with $\lambda`, present values of sums of moving-average coefficients are the same for returns and cash-flow growth.{^footnotechallenges} This present-value restriction {eq}`eqn:CSpresentval` in turn implies that

```{math}
:label: eqn:CampShill_punch
\alpha_v(\zeta)  = {\frac { \alpha_g(\zeta) -  \alpha_g(\lambda)}{\zeta - \lambda}}  -
{\frac { \alpha_r(\zeta) -  \alpha_r(\lambda)}{\zeta - \lambda}} .
```

The first term on the right side of {eq}`eqn:CampShill_punch` is the contribution to $\alpha_v(\zeta)$ of the expected cash flow in the absence of expected return variability, while the second  term accounts for expected  return variability.
In the time domain, equation {eq}`eqn:CampShill_punch` asserts  that

```{math}
\log V_t - \log G_t - \mu_v  = E \left[ \sum_{j=0}^\infty \lambda^j
\left( \log G_{t+1+j} - \log G_{t+j} - \mu_g - \log R_{t+1+j} + \mu_r \right) \vert {\mathfrak F}_t \right]
```

{cite:t}`CampbellShiller1988` use formulas like {eq}`eqn:CampShill_punch` to conclude that since

```{math}
E \left[ \sum_{j=0}^\infty \lambda^j
\left( \log G_{t+1+j} - \log G_{t+j} - \mu_g \right) \vert {\mathfrak F}_t \right]
```

has a small variance when $G_t$ is measured by aggregate dividends, observed large variations in the logarithm of the price-dividend ratio must be attributed to time-varying expected returns that act through the second term

```{math}
E \left[ \sum_{j=0}^\infty \lambda^j
\left( \log R_{t+j} - \mu_r \right) \vert {\mathfrak F}_t \right] .
```

That finding poses the scientific challenge of inventing a theoretical structure that generates variations in discounted expected returns that are not induced by dividend variations. We take that up in chapters XXXX.

[^footnotelabel]: See [Cointegration](sec:cointegrate) for a definition of cointegration.

[^footnotechallenges]: {cite:t}`HSRoberds` discuss challenges to testing this restriction. -->
