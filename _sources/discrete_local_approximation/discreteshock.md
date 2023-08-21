# Discrete Shock Elasticity

At the height of 

We use an approximation method to explore implications of the recursive
utility preference specification of [@KrepsPorteus:1978] and
[@EpsteinZin:1989] and counterparts to these preferences that capture
concerns about model misspecification. We present formulas for
(nonstandard) first and second-order approximations to dynamic,
stochastic equilibria for models in which economic agents have such
recursive preferences. The approximations build formulations from
[@schmittgrohe_uribe:2004] and [@LombardoUhlig:2018], we extend them in
a that features the uncertainty contributions more prominently. By
design, the implied approximations of stochastic discount factors used
to represent market or shadow values reside within the exponential
linear quadratic class. This class is known to give tractable formulas
for asset valuation over alternative investment horizons. See, for
instance, [@AngPiazzesi:2003] and [@Borovicka_Hansen:2013]. Moreover,
they are applicable to production-based macro-finance models with
investment opportunities in alternative forms of capital.

We use the approximations to provide further understanding of the
preferences and their implications for asset pricing. This opens the
door as well to other connections in the macroeconomics-finance
literature in which productive, investment and capital accumulation are
central model ingredients. As a central part of our analysis, we capture
the important uncertainty preference contribution as a change in the
probability distribution of the underlying economic dynamics. We link
this change of measure to the robust preferences specifications of
[@HansenSargent:2001] and [@AndersonHansenSargent:2003]. The robust
preferences formulations build on a robust control literature initiated
by [@Jacobson:1973] and [@Whittle:1981].

# Small noise expansion of the state dynamics {#sec:smallnoise101}

We follow [@LombardoUhlig:2018] by considering the following class of
stochastic processes indexed by a scalar perturbation parameter
$\mathsf{q}$:[^2]

$$X_{t+1}\left( \mathsf{q}\right) =\psi \left[ X_{t}\left( \mathsf{q}\right) ,%
\mathsf{q}W_{t+1},\mathsf{q}\right] . \label{eq:xlom}$$ Here $X$ is an
$n$-dimensional stochastic process and $\{W_{t+1}\}$ is an
i.i.d. normally distributed random vector with conditional mean vector
$0$ and conditional covariance matrix $I$. We parameterize this family
so that ${\sf q} = 1$ gives the model of interest.

We denote a zero-order expansion ${\sf q} = 0$ limit as:
$$\label{eqn:exporder0}
X^{0}_{t+1}  = \psi \left( X_{t}^0 ,0,0\right),$$ and assume that there
exists a second-order expansion of $X_{t}$ around $\mathsf{q}%
=0$: $$\label{eqn:exporder2}
X_{t} \approx X_{t}^0+\mathsf{q}X_{t}^1+\frac{{\sf q} ^2}{2}X_{t}^2$$
where $X_t^1$ is a first-order contribution and $X_t^2$ is a
second-order contribution.

In the remainder of this chapter we shall construct instances of the
second-order expansion
[\[eqn:exporder2\]](#eqn:exporder2){reference-type="eqref"
reference="eqn:exporder2"} in which the generic random variable $X_t$ is
replaced, for example, by the logarithm of consumption, a value
function, and so on. In approximation
[\[eqn:exporder2\]](#eqn:exporder2){reference-type="eqref"
reference="eqn:exporder2"}, the stochastic processes $X^j$, $j=0,1,2$
are appropriate derivatives of $%
X$ with respect to the perturbation parameter ${\sf q}$.

Processes $X_t^j, j=0, 1, 2$ have a recursive structure: the stochastic
process $X_{t}^0$ can be computed first, then the process $X_{t}^1$ next
(it depends on $X_t^0$), and then the process $X_{t}^2$ (it depends on
both $X_{t}^0$ and $X_{t}^1$).

We use a prime ($'$) to denote a transpose of a matrix or vector. When
we include $x'$ in a partial derivative of a scalar function it means
that the partial derivative is a row vector. Consistent with this
convention, let $\psi^i_{x'}$, the $i^{th}$ entry of $\psi_{x'}$, denote
the row vector of first derivatives with respect to the vector $x$, and
similarly for $\psi^i_{w'}$. Since ${\sf q}$ is scalar, $\psi_{\sf q}^i$
is the scalar derivative with respect to ${\sf q}$. Derivatives are
evaluated at $X_t^0$, which in many examples is invariant over time,
unless otherwise stated. This invariance follows when we impose a steady
state on the deterministic system.

The first-derivative process obeys a recursion $$\label{exprecur1}
X_{t+1}^1 =  \begin{bmatrix} \psi_{x'}^1 \cr \psi_{x'}^2 \cr \vdots \cr \psi_{x'} ^n \end{bmatrix} X_t^1 + \begin{bmatrix} \psi_{w'}^1 \cr \psi^2_{w'}
\cr \vdots \cr \psi^n_{w'} \end{bmatrix} W_{t+1} + \begin{bmatrix} \psi^1_{\sf q}  \cr \psi^2_{\sf q}  \cr \vdots \cr \psi^n_{\sf q}  \end{bmatrix}$$
that we can write compactly as the following *a first-order vector
autoregression*:
$$X_{t+1}^1 = \psi _{x'}X_{t}^1+\psi _{w'}W_{t+1}+\psi_{\sf q}$$ We
assume that the matrix $\psi_x'$ is stable in the sense that all of its
eigenvalues are strictly less than one in modulus.

It is natural for us to denote second derivative processes with double
subscripts. For instance, for the double script used in conjunction with
the second derivative matrix of $\psi^i$, the first subscript without a
prime ($'$) reports the row location; second subscript with a prime
($'$) reports the column location. Differentiating recursion
[\[exprecur1\]](#exprecur1){reference-type="eqref"
reference="exprecur1"} gives: $$\begin{aligned}
\label{exprecur2}
X_{t+1}^2 =  \hspace{.2cm} & \psi_{x'} X_t^2 + \begin{bmatrix} X_t^{1'} \psi_{xx'}^1 X_t^1 \cr X_t^{1'} \psi_{xx'}^2  X_t^1\cr \vdots \cr X_t^{1'} \psi_{xx'}^n X_t^1 \end{bmatrix} +  2 \begin{bmatrix} X_t^{1'} \psi_{xw'}^1 W_{t+1}  \cr X_t^{1'} \psi_{xw'}^2  W_{t+1} \cr \vdots \cr X_t^{1'} \psi_{xw'}^n W_{t+1}  \end{bmatrix} + \begin{bmatrix} {W_{t+1}}' \psi_{ww'}^1 W_{t+1}  \cr {W_{t+1}}' \psi_{ww'}^2  W_{t+1} \cr \vdots \cr {W_{t+1}}' \psi_{ww'}^n W_{t+1}  \end{bmatrix} \cr
& +  2 \begin{bmatrix} \psi^1_{{\sf q}  x'}X_{t}^1 \cr  \psi^2_{{\sf q}  x'}X_{t}^1 \cr \vdots \cr \psi^n_{{\sf q}  x'} X_t^1 \end{bmatrix}
+ 2 \begin{bmatrix} \psi_{{\sf q}  w' }^1W_{t+1} \cr \psi_{{\sf q}  w' }^2 W_{t+1} \cr \vdots \cr \psi^n_{{\sf q}  w'} W_{t+1} \end{bmatrix}
+ \begin{bmatrix} \psi^1_{{\sf q}  {\sf q} } \cr \psi^2_{{\sf q} {\sf q} } \cr \vdots \cr \psi^n_{{\sf q} {\sf q} } \end{bmatrix}
\end{aligned}$$ Recursions
[\[exprecur1\]](#exprecur1){reference-type="eqref"
reference="exprecur1"} and
[\[exprecur2\]](#exprecur2){reference-type="eqref"
reference="exprecur2"} have a linear structure with some notable
properties. The law of motion for $X^0$ is deterministic and is time
invariant if [\[eq:xlom\]](#eq:xlom){reference-type="eqref"
reference="eq:xlom"} comes from a stationary $\{X_t\}$ process. The
dynamics for $X^2$ are nonlinear only in $X^1$ and $W_{t+1}$. Thus, the
stable dynamics for $X^1$ that prevail when $\psi_x$ is a stable matrix
imply stable dynamics for $X^2$.

Let $C$ denote consumption and ${\widehat C}$ the logarithm of
consumption. Suppose that the logarithm of consumption evolves as:
$${\widehat C}_{t+1}  - {\widehat C}_t = \kappa (X_t, {\sf q} W_{t+1}, {\sf q}).$$
Approximate this process by: $$\label{eqn:exporder2c}
{\widehat C}_{t+1} - {\widehat C}_t \approx 
{\widehat C}_{t+1}^0 - {\widehat C}_t^0 + {\sf{q}} \left({\widehat C}_{t+1}^1 - {\widehat C}_t^1 \right)
+ \frac{{\sf q}^2} 2 \left({\widehat C}_{t+1}^2 - {\widehat C}_t^2\right)$$
where $$\begin{aligned}
{\widehat C}_{t+1}^0  - {\widehat C}_t^0  = \hspace{.2cm}  & \kappa (X_t^0, 0, 0) \eqdef \eta_0^c \cr 
{\widehat C}_{t+1}^1  - {\widehat C}_t^1  = \hspace{.2cm} & \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1} + \kappa_q \cr 
{\widehat C}_{t+1}^2  - {\widehat C}_t^2  =   \hspace{.2cm}  & \kappa_{x'} X_t^2 + {X_{t}^1}' \kappa_{x,x'} X_t^1 + 2
{X_{t}^1}' \kappa_{xw'} W_{t+1} + {W_{t+1}}'\kappa_{ww'}W_{t+1} \cr
&+ 2 \kappa_{q,x'} X_t^1 +  2 \kappa_{qw'}W_{t+1}  + \kappa_{qq} .
\end{aligned}$$ In models with endogenous investment and savings, the
consumption dynamics and some of the state dynamics will emerge as the
solution to a dynamic stochastic equilibrium model. We use the
approximating processes
[\[eqn:exporder2\]](#eqn:exporder2){reference-type="eqref"
reference="eqn:exporder2"} and
[\[eqn:exporder2c\]](#eqn:exporder2c){reference-type="eqref"
reference="eqn:exporder2c"} as inputs into the construction of an
approximating continuation value process and its risk-adjusted
counterpart for recursive utility preferences.

# Approximating a recursive utility value function {#sec:approxrecurutil}

In this section, we construct second-order expansions for components of
a continuation value process. This process along with its associated
stochastic discount factor process are important constituents of models.

The homogeneous of degree one representation of recursive utility is
$$\label{homog1a}
V_t = \left[ (1 - \beta) \left(C_t\right)^{1-\rho}  + \beta  \left( R_t \right)^{1-\rho} \right]^{\frac 1 {1-\rho}}$$
where $$\label{homog1b}
R_t = \left( {\mathbb E} \left[ \left( V_{t+1} \right)^{1-\gamma} \mid {\mathfrak A}_t \right] \right)^{\frac 1 {1-\gamma}} .$$
Notice that in equation [\[homog1a\]](#homog1a){reference-type="eqref"
reference="homog1a"}, $V_t$ is a homogeneous of degree one function of
$C_t$ and $R_t$. In equation
[\[homog1b\]](#homog1b){reference-type="eqref" reference="homog1b"},
$R_t$ is a homogeneous of degree one function of another function,
namely, $V_{t+1}$ as it varies over date $t+1$ information. In equation
[\[homog1a\]](#homog1a){reference-type="eqref" reference="homog1a"},
$0 < \beta < 1$ is a subjective discount factor and $\rho$ describes
attitudes toward intertemporal substitution. Formally, ${\frac 1 \rho}$
is the elasticity of intertemporal substitution. In equation
[\[homog1b\]](#homog1b){reference-type="eqref" reference="homog1b"},
$\gamma$ describes attitudes towards risk.

Continuation values are determined only up to an increasing
transformation. For computational and conceptual reasons, we find it
advantageous to work with the logarithm ${\widehat V}_t = \log V_t$. The
corresponding recursions for ${\widehat  V}_t$ expressed in terms of the
logarithm of consumption ${\widehat  C}_t$ are $$\label{eq:value_recur5}
{\widehat  V}_t = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t] + \beta \exp \left[(1-\rho) {\widehat  R}_t \right] \right]$$
where $$\label{eq:value_risk6}
{\widehat  R}_t = {\frac 1 {1 - \gamma}} \log {\mathbb E} \left(  \exp \left[ (1 - \gamma) {\widehat V}_{t+1} \right] \mid {\mathfrak A}_t \right).$$
The right side of recursion
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} is the logarithm of a constant elasticity
of substitution (CES) function of $\exp({\widehat  C}_t)$ and
$\exp({\widehat  R}_t)$.

::: remark
**Remark 1**. *The limit of ${\widehat R}_t$ as $\gamma$ approaches $1$
is ordinary expected logarithmic utility:
$$\lim_{\gamma \downarrow 1} \widehat   R_t = \lim_{\gamma \downarrow 1} {\frac { \log E \left( \exp\left[
(1-\gamma) {\widehat V}_{t+1}\right]
\vert {\mathfrak A}_t \right)}{1-\gamma}} = {\mathbb E}\left( \widehat  V_{t+1} \vert {\mathfrak A}_t \right) .$$*
:::

Our approach will be to construct small noise expansions for both
${\widehat  V}_t$ and ${\widehat  R}_t$ and then to assemble them
appropriately. Before doing so, we consider a reinterpretation of
[\[eq:value_risk6\]](#eq:value_risk6){reference-type="eqref"
reference="eq:value_risk6"}.

## Robustness to Model Misspecification {#sec:robust}

A reinterpretation of the utility recursion and the small-noise
expansion approach that we'll deploy comes from recognizing that when
$\gamma > 1$,
[\[eq:value_risk6\]](#eq:value_risk6){reference-type="eqref"
reference="eq:value_risk6"} emerges from an instance robust control
theory in which $\frac 1 {\gamma - 1}$ is a penalty parameter on entropy
relative to alternatives that constrains the alternative probability
models that a decision maker considers when evaluating consumption
processes. This interpretation originated in work by [@Jacobson:1973]
and [@Whittle:1981] that was extended and reformulated recursively by
[@Hansen1995].

Let the random variable $N_{t+1}  \ge 0$ satisfy
${\mathbb E} \left( N_{t+1} \mid {\mathfrak A}_t \right) = 1$ so that it
is a likelihood ratio. Think of replacing the expected continuation
value
${\mathbb E} \left( {\widehat  V}_{t+1} \mid {\mathfrak A}_t \right)$ by
the minimized value of the following problem:
$$\label{eqn:minproblem101}
\min_{N_{t+1} \ge 0, {\mathbb E}\left(N_{t+1} \vert {\mathfrak A}_t \right) = 1} {\mathbb E} \left( N_{t+1} {\widehat  V}_{t+1} \mid {\mathfrak A}_t \right) + \xi  {\mathbb E} \left( N_{t+1} \log N_{t+1}  \mid {\mathfrak A}_t \right)$$
where $\xi$ is a parameter that penalizes departures of $N_{t+1}$ from
unity as measured by relative entropy. Conditional relative entropy for
an altered conditional probability induced by applying change of measure
$N_{t+1}$ is
$${\mathbb E} \left( N_{t+1} \log N_{t+1}  \mid {\mathfrak A}_t \right) \ge 0$$
where, because the function $y \log y$ is convex, the inequality follows
from Jensen's inequality. Relative entropy is zero when $N_{t+1} = 1$.

::: remark
**Remark 2**. *To solve minimization problem
[\[eqn:minproblem101\]](#eqn:minproblem101){reference-type="ref"
reference="eqn:minproblem101"}, introduce a Lagrange multiplier, $\ell,$
on the conditional expectation constraint. The Lagrangian problem
separates across states, leading to the unconstrained problem:
$$\min_n n v + \xi n \log n + \ell(n - 1)$$ where $n$ is a potential
realization of $N_{t+1}$ and $v$ is a realization of $V_{t+1}.$ The
first-order conditions are: $$v + \xi + \xi \log n + \ell = 0.$$ The
solution is
$$n^* = \exp\left[ - \frac 1 \xi \left(v + \ell + \xi \right) \right],$$
and the minimizing objective
$$- \xi \exp\left[ - \frac 1 \xi \left(v + \ell + \xi\right) \right] - \ell .$$
To complete the solution, we solve for $\ell,$
$$\max_\ell - \xi {\mathbb E} \left( \exp\left[ - \frac 1 \xi \left(V_{t+1} + \ell  + \xi \right) \right] \mid {\mathfrak A}_t \right) - \ell$$
The first-order conditions are:
$${\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right) \exp\left[ - \left(\frac{ \ell + \xi}   \xi\right) \right] -1 = 0.$$
Thus the solution for $\ell$ is
$$\ell^* = \xi \log  {\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right) - 1$$
with a minimized objective given by
$$-  \xi \log  {\mathbb E} \left( \exp\left[ - \left(\frac 1 \xi \right) V_{t+1}  \right] \mid {\mathfrak A}_t \right).$$
The implied minimizer for $N_{t+1}$ is $$\label{Nstar}
N_{t+1}^*  = {\frac { \exp \left( - \frac 1 \xi {\widehat  V}_{t+1} \right) }{ {\mathbb E} \left[ \exp \left( - \frac 1 \xi {\widehat  V}_{t+1} \right) \mid
{\mathfrak A}_t \right]}} .$$*
:::

The minimizer of problem
[\[eqn:minproblem101\]](#eqn:minproblem101){reference-type="eqref"
reference="eqn:minproblem101"} given by
[\[Nstar\]](#Nstar){reference-type="eqref" reference="Nstar"} "tilts"
probabilities towards low continuation values, a version of what
[@Bucklew:2004] calls a stochastic version of Murphy's law. Notice that
the minimized objective satisfies
$$- \xi \log {\mathbb E} \left[ \exp \left( - {\frac 1 \xi} {\widehat  V}_{t+1} \right) \mid
{\mathfrak A}_t \right] = {\widehat R}_t 
%= {\frac 1 {1 - \gamma}} \log {\mathbb E} \left[ \exp \left( - {\frac 1 \xi} {\widehat  V}_{t+1} \right) \mid
%{\mathfrak A}_t \right]$$ where ${\widehat  R}_t$ was given previously
by equation [\[eq:value_risk6\]](#eq:value_risk6){reference-type="eqref"
reference="eq:value_risk6"} if we set $\xi = {\frac 1 {\gamma - 1}}$.
The random variable $N_{t+1}^*$ will play a central role in the
discussion that follows.

## Our expansion protocol

To approximate the recursive utility process, we deviate from common
practice in macroeconomics by letting the risk aversion or robust
parameter in preferences depend on ${\sf q}:$
$$\xi = {\sf q} \xi_o    \hspace{1 cm}  \gamma - 1 =  \frac {\gamma_o - 1} {\sf q}$$
The aversion to model misspecification or the aversion to risk moves
inversely with the parameter ${\sf q}$ when we embed the model of
interest within a parameterized family of models. In effect, the
variable ${\sf q}$ is doing double duty. Reducing ${\sf q} > 0$ limits
the overall exposure of the economy to the underlying shocks. This is
offset by letting the preferences include a greater aversion to
uncertainty. This choice of any expansion protocol has significant and
enlightening consequences for continuation value processes and for the
minimizing $N$ process used to alter expectations. It has antecedents in
the control theory literature, and it has the virtue that implied
uncertainty adjustments occur more prominently at lower-order terms in
the approximation.

### Order-zero

Write the order-zero expansion of
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} as $$\begin{aligned}
{\widehat  V}_t^0 & = {\frac 1 {1 - \rho}} \log \left[ (1 - \beta) \exp[(1-\rho) {\widehat  C}_t^0] + \beta \exp \left[(1-\rho) {\widehat  R}_t^0 \right] \right] \cr
{\widehat  R}_t^0 &  =  {\widehat  V}_{t+1}^0,
\end{aligned}$$ where the second equation follows from noting that
randomness vanishes in the limit as $\mathsf{q}$ approaches $0$.

For order zero, write the consumption growth-rate process as
$${\widehat C}_{t+1}^0 - {\widehat  C}_t^0  = \eta_c^0 .$$ The
order-zero approximation of
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} is:
$${\widehat  V}_t^0 - {\widehat  C}_t^0 = {\frac 1 {1 - \rho}}  \log \left[ (1 - \beta)  + \beta \exp \left[(1-\rho) \left( {\widehat  V}_{t+1}^0 - {\widehat C}_{t+1}^0 + \eta_c^0 \right)  \right] \right]$$
We guess that ${\widehat  V}_t^0 - {\widehat  C}_t^0 =\eta_{v-c}^0$ and
will have verified the guess if the following equation is satisfied
$$\exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right] = (1 - \beta) + \beta \exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right]\exp \left[ (1 - \rho) \eta_c^0 \right],$$
which implies $$\label{zero_formula}
\exp\left[(1-\rho) {\left( \eta_{v - c}^0 \right) }\right] = {\frac {1 - \beta} { 1 - \beta \exp \left[ (1 - \rho) \eta_c^0 \right]}} .$$
Equation [\[zero_formula\]](#zero_formula){reference-type="eqref"
reference="zero_formula"} determines $\eta_{v - c}^0$ as a function of
$\eta_c^0$ and the preference parameters $\rho, \beta$, but not the risk
aversion parameter $\gamma$. Specifically, $$\label{zero_formula2}
 \eta_{v-c}^0 = \frac {\log (1 - \beta)  - \log \left(1 - \beta \exp \left[ (1 - \rho) \eta_c^0 \right]\right)}{1 - \rho}$$

### Order-one

We temporarily take ${\widehat  R}_t^1 - {\widehat C}_t^1$ as given
(we'll compute it in section
[\[sec:second_order\]](#sec:second_order){reference-type="eqref"
reference="sec:second_order"}). We construct a recursion by taking a
first-order approximation to the nonlinear utility recursion
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} $$\label{eqn:lambda1}
{\widehat  V}_t^1 - {\widehat C}_t^1 = \lambda \left({\widehat  R}_t^1 - {\widehat C}_t^1\right)$$
where $$\begin{aligned}
\label{eqn:lambda2}
\lambda  & = \left[ {\frac {\beta \exp \left[(1-\rho) \left(\eta_{ v - c}^0 +\eta_c^0 \right ) \right]} { (1 - \beta)  + \beta \exp \left[(1-\rho) \left(\eta_{ v - c}^0 + \eta_c^0  \right ) \right]}} \right] \cr
& = \left[ {\frac {\beta \exp \left[(1-\rho) \eta_c^0\right]} { (1 - \beta) \exp\left[- (1-\rho) \eta_{ v - c}^0\right]   + \beta \exp \left[(1-\rho) \eta_c^0 \right]}} \right] \cr
& = \left[ {\frac {\beta \exp \left[(1-\rho) \eta_c^0 \right]} {  1 - \beta \exp \left[ (1 - \rho)\eta_c^0 \right]  + \beta \exp \left[(1-\rho) \eta_c^0 \right]}} \right] \cr
& = \beta \exp \left[(1-\rho) \eta_c^0 \right]
\end{aligned}$$ Notice how parameter $\rho$ influences the weight
$\lambda$ when $\eta_c \neq 0$, in which case the log consumption
process displays growth or decay. When $0 < \rho < 1$, the condition
$\lambda  <1$ restricts the parameter $\rho$ relative to the consumption
growth rate $\eta_c$ since $$(1- \rho) \eta_c <  - \log \beta$$

To facilitate computing some useful limits we construct:
$$\begin{aligned}
 \label{tilde_transform}
{\widetilde V}_t & = {\frac {{\widehat  V}_t - {{\widehat  V}_t^0}}{{\sf q} }} \cr
{\widetilde R}_t & =  {\frac {{\widehat  R}_t - {{\widehat  V}_{t+1}^0}}{{\sf q} }}
\end{aligned}$$ which we assume remain well defined as ${\sf q}$
declines to zero, with limits denoted by
${\widetilde V}_t^0, {\widetilde  R}_t^0$. Importantly,
$$\label{tilde_R}
{\widetilde R}_t =  \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right] \mid {\mathfrak A}_t \right),$$
Taking limits as ${\sf q}$ declines to zero:
$${\widehat R}_t^1 =  \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widehat V}_{t+1}^1 \right] \mid {\mathfrak A}_t \right)$$
Subtracting ${\widehat C}_t^1$ from both sides gives:
$$\label{first_order_risk}
{\widehat R}_t^1 - {\widehat C}_t^1 =  \left( {\frac 1 {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1} ^1 \right) + (1-\gamma_o)  \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t} ^1 \right)  \right] \mid {\mathfrak A}_t \right)$$
Substituting formula
[\[first_order_risk\]](#first_order_risk){reference-type="eqref"
reference="first_order_risk"} into the right side of
[\[eqn:lambda1\]](#eqn:lambda1){reference-type="eqref"
reference="eqn:lambda1"} gives the recursion for the first-order
continuation value: $$\label{first_recursive_update}
{\widehat V}_t^1 - {\widehat C}_t^1 =   \left( {\frac \lambda {1 - \gamma_o}} \right) \log {\mathbb E} \left( \exp \left[ (1 - \gamma_o)  \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1} ^1 \right) + (1-\gamma_o)  \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t} ^1 \right)  \right] \mid {\mathfrak A}_t \right)$$

::: remark
**Remark 3**. *We produce a solution by "guess and verify." Suppose that
$$\label{first-solution}
{\widehat V}_t^1 - {\widehat C}_t^1 = {\upsilon_1}' X_t^1 + \upsilon_0$$
It follows from
[\[first_recursive_update\]](#first_recursive_update){reference-type="eqref"
reference="first_recursive_update"} that $$\begin{aligned}
 \label{first_formula}
{\upsilon_1}' = & \lambda \left({\upsilon_1}' \psi_{x'}  +  \kappa_{x'} \right)\cr
\upsilon_0 = & \lambda \left( \upsilon_0 +   {\upsilon_1}'  \psi_{q} +  \kappa_q + {\frac {(1 - \gamma_0)} 2} 
\left| {\upsilon_1}' \psi_{w'}  
+ \kappa_{w'} \right|^2 \right).  
\end{aligned}$$ Deduce the second equation by observing that
$\exp \left[ (1 - \gamma_o)  \left( {\widehat V}_{t+1}^1 - {\widehat C}_{t+1} ^1 \right) + (1-\gamma_o)  \left( {\widehat C}_{t+1}^1 - {\widehat C}_{t} ^1 \right) \right]$
is distributed as a log normal. The solutions to equations
[\[first_formula\]](#first_formula){reference-type="eqref"
reference="first_formula"} are: $$\begin{aligned}
\upsilon_1 = & \lambda  \left( I - \lambda \psi_{x'} \right)^{-1} \kappa_{x'} \cr
\upsilon_0 = & {\frac \lambda {(1 - \lambda)}} \left({\upsilon_1}'  \psi_{q} +  \kappa_{q} \right) +  
 \frac {\lambda(1 - \gamma_0)} {2(1 - \lambda)} 
\left| {\upsilon_1}' \psi_{w'} + \kappa_{w'}  \right|^2.  
\end{aligned}$$ The continuation value has two components. The first is:
$${\upsilon_1}'X_t^1  + {\frac \lambda {(1 - \lambda)}} \left({\upsilon_1}'  \psi_{q} +  \kappa_{q} \right)
= {\mathbb E} \left[ \sum_{j=1}^\infty \lambda^j \left({\widehat C}_{t+j}^1 - {\widehat C}_{t+j-1} ^1\right) \mid {\mathfrak A}_t  \right]$$
and the second is a constant long-run risk adjustment given by:
$$\frac {\lambda(1 - \gamma_o)} {2(1 - \lambda)} 
\left| {\upsilon_1}' \psi_{w'} + \kappa_{w'} \right|^2.$$ This second
term is the the variance of $$\label{discount_future}
 {\mathbb E} \left[\sum_{j=1}^\infty \lambda^j \left({\widehat C}_{t+j }^1 - {\widehat C}_{t+j-1} ^1\right) \mid {\mathfrak A}_{t+1}  \right]$$
conditioned on ${\mathfrak A}_t$ scaled by
$\frac {\lambda(1 - \gamma_o)} {2(1 - \lambda)}$.*
:::

::: remark
**Remark 4**. *The formula for $\upsilon_1$ depends on the parameter
$\rho$. Moreover, $\upsilon_1$ has a well defined limit as $\lambda$
tends to unity as does the variance of
[\[discount_future\]](#discount_future){reference-type="eqref"
reference="discount_future"}. This limiting variance:
$$\lim_{\lambda \rightarrow 1} \left| {\upsilon_1}' \psi_{w'} + \kappa_{w'} \right|^2.$$
converges to the variance of the martingale increment of
${\widehat C}^1$.*
:::

::: {#remark:first .remark}
**Remark 5**. *Consider the logarithm of the risk adjusted continuation
value approximated to the first order. Note that from
[\[first-solution\]](#first-solution){reference-type="eqref"
reference="first-solution"},
$${\widehat V}_{t+1}^1 - {\widehat C}_t^1  = {\upsilon_1}' X_{t+1}^1 + \upsilon_0 + \kappa_{x'} X_t^1 + \kappa_{w'} W_{t+1} .$$
Substitute this expression into fromula
[\[first_order_risk\]](#first_order_risk){reference-type="eqref"
reference="first_order_risk"} and use the formula for the mean of random
variable distributed as a log normal to show that
$${\widehat V}_{t+1}^1- {\widehat R}_t^1 = \left( {\upsilon_1}'\psi{_w'} + \kappa_{w'} \right) W_{t+1} - \left( \frac { 1 - \gamma_o } 2 \right)
\vert  {\upsilon_1}'\psi{_w'} + \kappa_{w'} \vert^2$$*
:::

Equation
[\[first_order_risk\]](#first_order_risk){reference-type="eqref"
reference="first_order_risk"} is a standard risk-sensitive recursion
applied to log-linear dynamics. For instance, see [@Tallarini:2000]'s
paper on risk-sensitive business cycles and [@HansenHeatonLi:2008]'s
paper on measurement and inference challenges created by the presence of
long-term risk. Both of those papers assumed a logarithmic one-period
utility function, so that for them $\rho=1.$

Here we have instead obtained the recursion as a first-order
approximation without necessarily assuming log utility. Allowing for
$\rho$ to be different than one shows up in both the order zero and
order one approximations as reflected in
[\[zero_formula2\]](#zero_formula2){reference-type="eqref"
reference="zero_formula2"} and
[\[first_recursive_update\]](#first_recursive_update){reference-type="eqref"
reference="first_recursive_update"}, respectively. As reflected by
formula
[\[first_recursive_update\]](#first_recursive_update){reference-type="eqref"
reference="first_recursive_update"}, for the first-order approximation
the parameter $\lambda = \beta$ when $\rho = 1$. But otherwise, it is
different. Equation
[\[first_order_risk\]](#first_order_risk){reference-type="eqref"
reference="first_order_risk"} also is very similar to a first-order
approximation proposed in [@RestoyWeil:2011]. Like formula
[\[first_order_risk\]](#first_order_risk){reference-type="eqref"
reference="first_order_risk"}, @RestoyWeil:2011 allow for $\rho \ne 1$.
In contrast, our equation has an explicit constant term coming from the
risk/robustness adjustment, and we have explicit formula for $\lambda$
that depends on preference parameters and the consumption growth rate.

### Order two {#sec:second_order}

Differentiating equation
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} a second time gives: $$\begin{aligned}
{\widehat  V}_t^2  = 
%& \hspace{.2cm} (1 - \lambda)  {\widehat  C}_t^2 +  \lambda  {\widehat  R}_t^2 \cr & +
%(1- \rho)\left[ (1 - \lambda)\left({\widehat  C}_t^1\right)^2 +  \lambda  \left( {\widehat  R}_t^1 \right)^2 - \left[ (1 - \lambda) {\widehat  C}_t^1  + \lambda  {\widehat  R}_t^1 \right]^2 \right] \cr
 & \hspace{.2cm} (1 - \lambda)  {\widehat  C}_t^2 +  \lambda  {\widehat  R}_t^2  + (1- \rho) (1-\lambda) \lambda  \left( {\widehat  R}_t^1 - {\widehat  C}_t^1\right)^2 .
\end{aligned}$$ Equivalently,
$${\widehat  V}_t^2 - {\widehat  C}_t^2 = \lambda \left( {\widehat  R}_t^2 - {\widehat  C}_t^2 \right) +  (1- \rho) (1-\lambda) \lambda  \left( {\widehat  R}_t^1 - {\widehat  C}_t^1\right)^2 .$$

Rewrite transformation
[\[tilde_transform\]](#tilde_transform){reference-type="eqref"
reference="tilde_transform"} as $$\begin{aligned}
{\sf q} {\widetilde V}_t & = {\widehat  V}_t - {{\widehat  V}_t^0} \cr
{\sf q} {\widetilde R}_t & =  {\widehat  R}_t - {{\widehat  V}_{t+1}^0}
\end{aligned}$$ Differentiating twice with respect to ${\sf q}$ and
evaluated at ${\sf q} = 0$ $$\begin{aligned}
\left. 2 \frac {d} {d{\sf q}}{\widetilde V}_t  + {\sf q}  \frac {d^2} {d{\sf q}^2} {\widetilde V}_t \right\vert_{{\sf q} = 0} & =  2{\widetilde  V}_t^1   = {\widehat  V}_t^2  \cr
\left. 2 \frac {d} {d{\sf q}}{\widetilde R}_t  + {\sf q}  \frac {d^2} {d{\sf q}^2} {\widetilde R}_t \right\vert_{{\sf q} = 0} & =  2{\widetilde  R}_t^1 =  {\widehat  R}_t^2
\end{aligned}$$ Differentiating
[\[tilde_R\]](#tilde_R){reference-type="eqref" reference="tilde_R"} with
respect to ${\sf q}$
$${\frac {d{\widetilde R}_t} {d {\sf q} }} =   {\frac {{\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right] {\frac {d{\widetilde V}_{t+1}} {d {\sf q} }} \mid {\mathfrak A}_t \right)}{{\mathbb E} \left( \exp \left[ (1 - \gamma_o)  {\widetilde V}_{t+1} \right]  \mid {\mathfrak A}_t \right)}},$$
and thus $$\begin{aligned}
 \label{recur_update}
{\widehat  R}_t^2 = 2 {\widetilde R}_t^1 & = 2 E \left( N_{t+1}^0  {\widetilde V}^1_{t+1} \mid {\mathfrak A}_t \right)  \cr
& =  E \left( N_{t+1}^0  {\widehat V}^2_{t+1} \mid {\mathfrak A}_t \right),
\end{aligned}$$ where $N_{t+1}^0$ $$\begin{aligned}
 \label{firstN}
N_{t+1}^0  \eqdef &  {\frac { \exp \left[ (1 - \gamma_o)   {\widetilde V}_{t+1}^0 \right]}{{\mathbb E}\left( \exp \left[ (1 - \gamma_o)   {\widetilde V}_{t+1}^0 \right] \mid {\mathfrak A}_t \right)}} \cr
= & {\frac { \exp \left[ (1 - \gamma_o)   {\widehat  V}_{t+1}^1 \right]}{{\mathbb E} \left( \exp \left[ (1 - \gamma_o)   {\widehat  V}_{t+1}^1 \right] \mid {\mathfrak A}_t \right)}} .
\end{aligned}$$ Subtracting ${\widehat C}_t^2$ from ${\widehat R}_t^2$
and substituting into
[\[recur_update\]](#recur_update){reference-type="eqref"
reference="recur_update"} gives: $$\label{second-order}
{\widehat  V}_t^2 - {\widehat  C}_t^2 = \lambda {\mathbb E} \left( N_{t+1}^0  \left[\left({\widehat V}^2_{t+1} - {\widehat C}^2_{t+1} \right)
+ \left( {\widehat C}^2_{t+1} -  {\widehat C}^2_{t}\right) \right] \mid {\mathfrak A}_t \right)
+  (1- \rho) (1-\lambda) \lambda  \left( {\widehat  R}_t^1 - {\widehat  C}_t^1\right)^2 .$$
Even if the second-order contribution to the consumption process is
zero, there will be nontrivial adjustment to the approximation of
${\widehat V} - {\widehat C}$ because
$\left( {\widehat  R}^1 - {\widehat  C}^1\right)^2$ is different from
zero. This term vanishes when $\rho = 1$, and its sign will be different
depending on whether $\rho$ is bigger or smaller than one.

::: remark
**Remark 6**. *The calculation reported in Remark
[5](#remark:first){reference-type="ref" reference="remark:first"}
implies that
$$\log N_{t+1}^0 = (1-\gamma_o)\left({\widehat V}_{t+1}^1- {\widehat R}_t^1\right) =   (1- \gamma_o) \left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right) W_{t+1} -  \frac { (1 - \gamma_o)^2 } 2 
\vert  {\upsilon_1}'\psi{_w'} + \kappa_{w'} \vert^2$$ As a consequence,
under the change in probability measure induced by $N_{t+1}^0,$
$W_{t+1}$ has a mean given by
$$\mu^0 \eqdef (1 - \gamma_o)  \left( {\upsilon_1}'\psi_{w'} + \kappa_{w'} \right)'$$
and with the same covariance matrix given by the identity. This is an
approximation to robustness adjustment expressed as an altered
distribution of the underlying shocks. It depends on
$\gamma_o -1 = {\frac 1 \xi_o}$ as well as the state dynamics as
reflected by $\upsilon_1$ and by the shock exposure vectors $\psi_{w'}$
and $\kappa_{w'}$.*
:::

# Stochastic discount factor process

A stochastic discount factor (SDF) process $S = \{ S_t : t \ge 0  \}$
tells how a consumer responds to small changes in uncertainty and
thereby consequently how a consumer values risky payouts. SDF processes
have a variety of uses. First, they provide shadow prices that tell how
a consumer's uncertainty aversion shapes marginal valuations of risky
assets. Second, they shape first-order conditions for optimally choosing
financial and physical investments. Third, they underly tractable
formulas for equilibrium asset prices. Fourth, they can help construct
Pigouvian taxes for correcting adverse externalities under uncertainty.
Fifth, they provide useful tools for assessing effects of small (local)
changes in government policies.

To indicate how to deduce an SDF process, we begin by positing that the
date zero value of a risky date $t$ consumption payout $\chi_t$ is
$$\label{eqn:price101}
\pi_0^t(\chi_t) = E\left[ \left( {\frac {S_t}{S_0}} \right) \chi_t  \Bigr| {\mathfrak A}_0 \right].$$
We compute the ratio $\frac {S_t}{S_0}$ that appears in formula
[\[eqn:price101\]](#eqn:price101){reference-type="eqref"
reference="eqn:price101"} by evaluating the slope of an indifference
curve that runs through both a baseline consumption process
$\{C_t\}_{t=0}^\infty$ and a perturbed consumption process
$$\left(C_0 - P_0({\sf q} ) , C_1, C_2, \ldots , C_t + {\sf q}  \chi_t, C_{t+1}, ... \right) .$$
We think of ${\sf q}$ as parameterizing an indifference curve, so
$P_0({\sf q} )$ expresses how much current period consumption must be
reduced to keep a consumer on the same indifference curve after we
replace $C_t$ by $C_t + {\sf q}  \chi_t$. We set $\pi_0^t(\chi_t)$
defined in equation
[\[eqn:price101\]](#eqn:price101){reference-type="eqref"
reference="eqn:price101"} equal to the slope of that indifference curve:
$$\pi_0^t(\chi_t) =  \left. {\frac d {d {\sf q}  } }  P_0({\sf q} ) \right|_{{\sf q}  = 0} .$$

The one-period increment in the stochastic discount factor process for
recursive utility is: $$\begin{aligned}
 \label{eqn:sdf50}
\frac {S_{t+1}}{S_t} & = \beta \left( \frac {C_{t+1} }{C_t} \right)^{-\rho}  \exp\left[ (1-\gamma) \left( {\widehat V}_{t+1} - {\widehat R}_t\right)\right] \exp\left[ (\rho- 1) \left( {\widehat V}_{t+1} - {\widehat R}_t\right)\right]\cr
&= \beta N_{t+1}^*   \exp \left( {\widehat S}_{t+1} - {\widehat S}_t \right)
\end{aligned}$$ where $${\widehat S}_{t+1} - {\widehat S}_t 
\eqdef \log \beta
- \rho\left( {\widehat C}_{t+1} + {\widehat C}_t\right)  + (\rho- 1) \left( {\widehat V}_{t+1} - {\widehat R}_{t} \right)$$
where $N_{t+1}^*$ induces the change of probability measure that we
described previously as the outcome of robustness problem. (See equation
[\[Nstar\]](#Nstar){reference-type="eqref" reference="Nstar"}.) We will
use this second formula in what follows.

::: remark
**Remark 7**. *To verify formula
[\[eqn:sdf50\]](#eqn:sdf50){reference-type="eqref"
reference="eqn:sdf50"}, we compute a one-period intertemporal marginal
rate of substitution. Given the valuation recursions
[\[eq:value_recur5\]](#eq:value_recur5){reference-type="eqref"
reference="eq:value_recur5"} and
[\[eq:value_risk6\]](#eq:value_risk6){reference-type="eqref"
reference="eq:value_risk6"}, we construct two marginal utilities
familiar from CES and exponential utility: $$\begin{aligned}
%mc & = (1 - \beta)\left( \frac  c  v \right)^{-\rho}  \cr
mc & = (1 - \beta)\left(   c   \right)^{-\rho} \exp\left[ (\rho - 1) {\hat v} \right]  \cr
%mr & = \beta \left( \frac  r v \right)^{1-\rho} \cr
m{\hat r} & = \beta  \exp[(1-\rho) \left({\hat r} - {\hat v}   \right) ]
\end{aligned}$$ From the certainty equivalent formula, we construct the
marginal utility of the next-period logarithm of the continuation value:
$$m{\hat v}^+=  \exp \left[(1-\gamma ) \left({\hat v}^+ -  {\hat r} \right) \right]$$
where the $+$ superscript is used to denote the next-period counterpart.
In addition, the next-period marginal utility of consumption is
$$mc^+ = (1 - \beta)\left(   c ^+  \right)^{-\rho} \exp\left[ (\rho - 1) {\hat v}^+\right]$$
Putting these four formulas together using the chain rule for
differentiation gives a marginal rate of substitution:
$$\frac {(mr) (mv^+) (mc^+)}{mc} = \beta \left( \frac {c^+}{c} \right)^{-\rho} \exp \left[(1-\gamma ) \left({\hat v}^+ -  {\hat r} \right) \right]
\exp \left[(\rho - 1) \left({\hat v}^+ -  {\hat r} \right) \right].$$
Now let ${\hat v}^+ =  {\widehat V}_{t+1}$, $c^+ = C_{t+1}$, $C_t = c$
and ${\hat r} = {\widehat R}_t$ to obtain the formula for the one-period
stochastic discount factor
[\[eqn:sdf50\]](#eqn:sdf50){reference-type="eqref"
reference="eqn:sdf50"}.*
:::

We approximate $\left[{\widehat S}_{t+1} - {\widehat S}_{t} \right]$ as
$${\widehat S}_{t+1} - {\widehat S}_{t}  \approx \left[{\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right] + \left[{\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right] + {\frac 1 2} 
\left[{\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right]$$ where
$$\begin{aligned}
{\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 & \eqdef  \log \beta - \rho \eta_c^0 \cr 
{\widehat S}_{t+1}^1 - {\widehat S}_{t}^1 & \eqdef - {\widehat C}_{t+1}^1 + {\widehat C}_{t}^1 + (\rho - 1) \left( {\widehat V}_{t+1}^1 - {\widehat C}^1_{t+1} \right)  - 
(\rho - 1) \left( {\widehat R}_t^1 - {\widehat C}_t^1  \right)
\cr
{\widehat S}_{t+1}^2 - {\widehat S}_{t}^2 & \eqdef - {\widehat C}_{t+1}^2 + {\widehat C}_{t}^2 + (\rho - 1) \left( {\widehat V}_{t+1}^2 - {\widehat C}^2_{t+1} \right)  - 
(\rho - 1) \left( {\widehat R}_t^2 - {\widehat C}_t^2  \right)
\end{aligned}$$

We now consider three different approaches to approximating $N_{t+1}^*$.

## Approach 1

Write $$\begin{aligned}
N_{t+1}^* & = \frac { \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right]}{{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right] \mid {\mathfrak A}_t \right]} \cr
& = \frac { \exp\left[(1-\gamma_o) {\widetilde V}_{t+1} \right]}{ \exp\left[(1-\gamma_o) {\widetilde R}_t\right] }  
\end{aligned}$$ Form the "first-order" approximation: $$\begin{aligned}
 \label{firsttilde}
\log N_{t+1}^*  & \approx (1 - \gamma_o) \left[\left( {\widetilde V}_{t+1}^0 - {\widetilde R}_t^0  \right) +   {\sf q} \left( {\widetilde V}_{t+1}^1 - {\widetilde R}_t^1 \right)\right] \cr
& = (1 - \gamma_o)\left[ \left( {\widehat V}_{t+1}^1 - {\widehat R}_t^1 \right) +  {\frac {\sf q}  2}  \left( {\widehat V}_{t+1}^2 - {\widehat R}_t^2 \right) \right]
\end{aligned}$$ This approach suggests using the following first-oder
approximation for the stochastic discount factor: $$\begin{aligned}
\log S_{t+1} -  &\log S_t \approx (1 - \gamma_o)\left[ \left( {\widehat V}_{t+1}^1 - {\widehat R}_t^1 \right) +  {\frac 1 2}  \left( {\widehat V}_{t+1}^2 - {\widehat R}_t^2 \right) \right] \cr & + \left[{\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right] + \left[{\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right]  
\end{aligned}$$ While the implied $N_{t+1}^*$ approximation is positive,
it will not have conditional expectation equal to one. In contrast, the
exponential of the first-order contribution
$(1 - \gamma_o) \left( {\widehat V}_{t+1}^1 -  {\widehat R}_t^1\right)$
will have conditional expectation equal to one as we have noted
previously.

## Approach 2

If we were to use a second-order approximation of $N_{t+1}^*$, it would
push us outside the class of exponentially quadratic stochastic discount
factors. Instead we could combine a first-order approximation of
$\log N_{t+1}^*$ with a second-order approximation of
${\widehat S}_{t+1} - {\widehat S}_t$: $$\begin{aligned}
\log S_{t+1} -  &\log S_t \approx (1 - \gamma_o)\left[ \left( {\widehat V}_{t+1}^1 - {\widehat R}_t^1 \right) +  {\frac 1 2}  \left( {\widehat V}_{t+1}^2 - {\widehat R}_t^2 \right) \right] \cr & + \left({\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right) + \left({\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right)  + \frac 1 2 \left({\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right) 
\end{aligned}$$ which would would preserve the quadratic approximation
of $\log S_{t+1} -  \log S_t.$

## Approach 3

Next consider an alternative modification of Approach 1 whereby:
$$\begin{aligned}
\log N_{t+1}^* & \approx \frac { \exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\widetilde V}_{t+1}^1 \right)\right]} 
{{\mathbb E}\left(\exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\widetilde V}_{t+1}^1 \right)\right] \mid {\mathfrak A}_t \right)} \cr
& = \frac { \exp\left[(1-\gamma_o)\left[ \left({\widehat V}_{t+1}^1 - {\widehat R}_t^1\right) + {\frac 1 2} \left( {\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right) \right]\right]} 
{{\mathbb E}\left(
\exp\left[(1-\gamma_o)\left[ \left({\widehat V}_{t+1}^1 - {\widehat R}_t^1\right) + {\frac 1 2} \left( {\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right) \right]\right]
\mid {\mathfrak A}_t \right)} . 
\end{aligned}$$ is used in conjunction with
$$\left({\widehat S}_{t+1}^0 - {\widehat S}_{t}^0 \right) + \left({\widehat S}_{t+1}^1 - {\widehat S}_{t}^1\right)  + \frac 1 2 \left({\widehat S}_{t+1}^2 - {\widehat S}_{t}^2\right) .$$
By design, exponential counterpart of this approximation will have
conditional expectation equal to one. With a little bit of algebraic
manipulation, it may be shown that this approximation induces a
distributional change for $W_{t+1}$ with a conditional mean that is
affine in $X_{t+1}$ and an altered conditional variance that is constant
over time.

To understand better this choice of approximation, consider the family
of random variables (indexed by ${\sf q}$) $$\label{family}
(1-\gamma_o) \left({\widetilde V}_{t+1}^0 + {\sf q} {\widetilde V}_{t+1}^1\right) 
- \log {\mathbb E} \left( \exp\left[(1-\gamma_o) \left({\widetilde V}_{t+1}^0 +  {\sf q} {\widetilde V}_{t+1}^1\right) \right] \mid {\mathfrak A}_t \right) .$$
The corresponding family of exponentials has conditional expectation one
and the ${\sf q} = 1$ member is the proposed approximation for
$N_{t+1}^*.$ Differentiate the family with respect to ${\sf q}$:
$${\widetilde V}_{t+1}^1 - \frac {{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1}^0  \right] {\widetilde V}_{t+1}^1 \mid {\mathfrak A}_t \right)}
{{\mathbb E} \left( \exp\left[(1-\gamma_o) {\widetilde V}_{t+1}^0  \right]  \mid {\mathfrak A}_t \right)} 
 =
{\widetilde V}_{t+1}^1 - {\widetilde R}_{t}^1.$$ Thus this family of
random variables has the same first-order approximation in ${\sf q}$ as
$\log N_{t+1}^*$ given in
[\[firsttilde\]](#firsttilde){reference-type="eqref"
reference="firsttilde"} and it remains within the linear-quadratic in
logarithms formulation.

As a first-order change of probability measure, this approximation will
induce state dependence in the conditional mean and will alter the
covariance matrix of the shock vector. We find this approach interesting
because it links back directly to the outcome of the robustness
formulation we described in Section
[3.1](#sec:robust){reference-type="ref" reference="sec:robust"}.
Moreover, the state dependence in the mean will induce a corresponding
state dependence in the one-period uncertainty prices.

# Long-run risk example

We consider a model with long-run risk components to consumption as
suggested by [@bansalyaron2004]. For the moment, we abstract from
production; but as we will see later there is a production counterpart
in consumption displays long-run risk. For now think simply specify a
consumption process with a long-run risk component.

## Approximation

By applying this approximation to the [@bansalyaron2004] model, we
obtain the state dynamics: $$\begin{aligned}
X_{t+1}^0 = &\begin{bmatrix} 0 \cr 1 \end{bmatrix} \cr
X_{t+1}^1 = & \begin{bmatrix}  \theta_{11}^x & 0 \cr 0 & \theta_{22}^x \end{bmatrix} X_t^1 + \begin{bmatrix} \sigma_{11}^x & 0 & 0 \cr 0 & \sigma_{22}^x & 0 \end{bmatrix} W_{t+1} \cr
X_{t+1}^2 = &\begin{bmatrix}  \theta_{11}^x & 0 \cr 0 & \theta_{22}^x \end{bmatrix}  X_t^2 +  \begin{bmatrix}  2X_{2,t}^1\sigma_{11}^x & 0 & 0 \cr 
0 & 0 & 0 \end{bmatrix} W_{t+1},
\end{aligned}$$ where $0 < \theta^{x}_{11} < 1$ and
$0 < \theta_{22}^x < 1$, and the consumption dynamics: $$\begin{aligned}
{\widehat C}_{t+1}^0 -  {\widehat C}_{t}^0 & = \eta_c^0\cr
{\widehat C}_{t+1}^1 - {\widehat C}_{t}^1 & =  \theta^c_1 X_{1,t}^1 + \begin{bmatrix} 0 & 0 & \sigma^c_3 \end{bmatrix} W_{t+1}\cr
{\widehat C}_{t+1}^2 - {\widehat C}_{t}^2 &  =  \theta^c_1 X_{1, t}^2 + \begin{bmatrix}  0 & 0 & 2X_{2,t}^1  \sigma^c_3 \end{bmatrix} W_{t+1}.  
\end{aligned}$$ Thus we take consumption to evolve (apprioximately) as:
$$\begin{aligned}
{\widehat C}_{t+1} - {\widehat C}_t & = \left({\widehat C}_{t+1}^0 + {\widehat C}_{t+1}^1 + {\frac 1 2} {\widehat C}_{t+1}^2\right)  -  \left({\widehat C}_{t}^0 +  {\widehat C}_{t}^1 + {\frac 1 2} {\widehat C}_{t}^2 \right) \cr
& = \eta_c^0 + \theta^c_1 \left(  X_{1,t}^1 +  {\frac 1 2} X_{1, t}^2 \right) + \begin{bmatrix} 0 & 0 & \left(1 +  X_{2,t}^1 \right) \sigma^c_3 \end{bmatrix} W_{t+1} . 
\end{aligned}$$ The processes $\{ X^1_{1,t} \}$ and $\{ X^2_{1,t} \}$
contribute temporal dependence to the consumption growth dynamics. The
process $\{ X_{2,t}^1 \}$ contributes stochastic volatility to the
consumption dynamics while the stationary specification of the process
$\{ X^2_{2,t} \}$ is identically zero and can be ignored.

## VAR approach

The Markov process governing the predictable component of macroeconomic
growth is scalar in the [@bansalyaron2004] analysis. Motivated by
empirical evidence, [@HansenHeatonLi:2008] study an extension of this
model where $X^1$ is a vector autoregression. To relate to the VAR
approach of @HansenHeatonLi:2008, write the first-order approximation
for the logarithm of consumption as:
$$\widehat C_{t+1}^1 - \widehat  C_t = \eta_c^0  + {\mathbb D} X_t^1 + {\mathbb F}' W_{t+1}$$
where $$\begin{aligned}
\kappa_{x'} &= {\mathbb D} \cr
\kappa_{w'} & = {\mathbb F}.
\end{aligned}$$ Where the first-order process $X^1$ includes a
predictable component of the macroeconomic growth-rate process and
evolves as an autoregression:
$$X_{t+1}^1 = {\mathbb A} X_t^1 + {\mathbb B} W_{t+1},$$ where
$$\begin{aligned}
\psi_{x'} & = {\mathbb A} \cr
\kappa_{w'} & = {\mathbb B}. 
\end{aligned}$$ and ${\mathbb A}$ is a stable matrix. Thus the
first-order approximation to the [@bansalyaron2004] for the consumption
dynamics is a special case of the formulation in
[@HansenHeatonLi:2008].[^3]

The row vector ${\mathbb F}$ and matrix ${\mathbb B}$ are configured so
that the components of the shock vector, $W_{t+1},$ directly disturbs
growth in the logarithm of consumption and its predictable (first-order)
growth component $X^1$. Notice, in particular that the conditional mean
of $\widehat  C_{t+j} - \widehat C_t$ is
$$j \eta_c^0 +{\mathbb D}\left( X_t + {\mathbb A} X_t + ... + {\mathbb A}^{j-1}\right) X_t.$$
The corresponding multi-period forecast errors contribute to the
variance of $\widehat  C_{t+j} - \widehat C_t$ with a variance that
increases with the horizon. When the process $\{ {\mathbb D} X_t\}$ is
highly persistent, there is said to be substantial "long-run risk" in
consumption.

## Approximating continuation values

Returning the original [@bansalyaron2004] specification, we consider the
approximation of continuation values and the corresponding change in
probabilities. The first-order continuation-value approximation is
$$\begin{aligned}
{\widehat V}_t^1 - {\widehat C}_t^1  &  = \left({\frac {\lambda}{1 - \lambda \theta_{11}^x}}  \right)\theta_1^c X_{1,t}^1 \cr
{\widehat R}_t^1 - {\widehat C}_t^1  &  = \left({\frac {1}{1 - \lambda \theta_{11}^x}}  \right)\theta_1^c X_{1,t}^1 \
\end{aligned}$$ The implied change in probability measure is
$$N_{t+1}^0  =  \exp\left(   \mu^0 \cdot W_{t+1}  - {\frac 1 2}  \left\vert \mu^0 \right\vert^2  \right)$$
where
$$\mu^0 = (1-\gamma_o) \begin{bmatrix} \left({\frac {\lambda}{1 - \lambda \theta_{11}^x}}  \right)\theta^c_1 \sigma_{11}^x \cr
0 \cr \sigma^c_3 \end{bmatrix}$$ is the implied mean distortion. The
negative of $\mu_0$ gives the vector of one period shock exposure
prices.

We use formula, [\[second-order\]](#second-order){reference-type="eqref"
reference="second-order"}, for the second-order adjustment. As a first
step we compute $${\mathbb E} \left[ N_{t+1}^0  
 \left( {\widehat C}^2_{t+1} -  {\widehat C}^2_{t}\right) \mid {\mathfrak A}_t \right] = \theta_c X_{1, t}^2 + 2 (\sigma_3)^2 X_{2,t}^1.$$
Thus we are lead to solve: $$\begin{aligned}
{\widehat  V}_t^2 - {\widehat  C}_t^2 = & \lambda {\mathbb E} \left[ N_{t+1}^1  \left({\widehat V}^2_{t+1} - {\widehat C}^2_{t+1} \right)
 \mid {\mathfrak A}_t \right] \cr & + \theta_c X_{1, t}^2 + 2 (\sigma_3)^2 X_{2,t}^1 + (1- \rho) (1-\lambda) \lambda\left[\left({\frac {1}{1 - \lambda \theta_{11}^x}}  \right)\theta_1^c X_{1,t}^1\right]^2 
\end{aligned}$$ forward using the change of probability measure under
which the conditional expectation of $W_{t+1}$ is equal to $\mu^0$.

## Shock elasticities

We use the shock elasticities to explore pricing implications of this
recursive utility specification. We conduct this exploration using the
original parameter calibration in [@bansalyaron2004]. These elasticities
and their relation to impulse-response functions introduced first to
macroeconomics by [@frisch:1933] is described
[@BorovickaHansenScheinkman:2014]. In what follows, we use
exponential/linear/quadratic implementation by [@Borovicka_Hansen:2013]
and by [@BorovickaHansen:2016].

Figure [1](#fig:expose){reference-type="ref" reference="fig:expose"}
gives the shock exposure elasticities for consumption to each of the
three shocks. This can can interpreted as nonlinear local impulse
responses for consumption (in levels not logarithms). The elasticities
for the growth rate shock and the stochastic volatility shock start
small and increase over the time horizon as dictated by the persistence
of the two exogenous state variable processes. The elasticities for the
direct shock to consumption are flat over the horizon as to be expected
since the shock directly impacts log consumption in a manor that is
permanent. Notice that while elasticities for the volatility shock are
different from zero, their contribution is much smaller than the other
shocks.[^4] Nevertheless, for this [@bansalyaron2004] calibration of the
long run risk model, stochastic volatility induces state dependence in
the elasticities for growth rate and consumption shocks as reflected by
quantiles given in the figures.

<figure id="fig:expose">

<figcaption>Exposure elasticities for three shocks. The time scale is in
months. </figcaption>
</figure>

Figure [2](#fig:price1){reference-type="ref" reference="fig:price1"}
gives the corresponding shock price elasticities for $\rho = 2/3$ and
$\gamma = 10$. The recursive utility preferences are forward looking as
reflected by the continuation-value contribution to the one-period
increment to the stochastic discount factor process as given in
[\[eqn:sdf50\]](#eqn:sdf50){reference-type="eqref"
reference="eqn:sdf50"}. This forward-looking contribution is reflected
in shock price elasticities that are now flat for both the growth rate
shock and the shock to stochastic volatility. The magnitudes are
substantially higher for the shock-price elasticities. While the
relative magnitudes are very different, the shock price elasticities are
much smaller than the other elasticities.[^5].

<figure id="fig:price1">

<figcaption>Price elasticities for three shocks. <span
class="math inline"><em>ρ</em> = 2/3, <em>γ</em> = 10,</span> <span
class="math inline"><em>β</em> = .998.</span> The time scale is in
months. </figcaption>
</figure>

Figures [3](#fig:price2){reference-type="ref" reference="fig:price2"}
and [4](#fig:price3){reference-type="ref" reference="fig:price3"}
provide the analogous plots for $\rho = 1, 1.5.$ The shock price
elasticities are very similar given these modes increases in $\rho$. It
is evidently the risk aversion parameter $\gamma = 10$ that is important
for determining the magnitude of these elasticities. Figure
[5](#fig:price4){reference-type="ref" reference="fig:price4"} sets
$\rho=\gamma = 10$ which corresponds to preferences that are time
separable. The forward-looking component to the stochastic discount
factor is shut down as is evident from formula
[\[eqn:sdf50\]](#eqn:sdf50){reference-type="eqref"
reference="eqn:sdf50"}. Now the shock price elasticities and shock
exposure elasticities show a very similar trajectory except that the
shock price elasticities are about ten times larger. The stochastic
volatility shock price is increased by about seventy-five times. Notice
that for longer horizons the $\gamma = \rho = 10$ preference model has
prices that are very similar in magnitude to the recursive preference
models with more modest specifications of $\rho$.

<figure id="fig:price2">

<figcaption>Price elasticities for three shocks. <span
class="math inline"><em>ρ</em> = 1, <em>γ</em> = 10,</span> <span
class="math inline"><em>β</em> = .998.</span> The time scale is in
months. </figcaption>
</figure>

<figure id="fig:price3">

<figcaption>Price elasticities for three shocks. <span
class="math inline"><em>ρ</em> = 3/2, <em>γ</em> = 10,</span> <span
class="math inline"><em>β</em> = .998.</span> The time scale is in
months. </figcaption>
</figure>

<figure id="fig:price4">

<figcaption>Price elasticities for three shocks. <span
class="math inline"><em>ρ</em> = 10, <em>γ</em> = 10,</span> <span
class="math inline"><em>β</em> = .998.</span> The time scale is in
months. </figcaption>
</figure>

# Solving models

The [@bansalyaron2004] example along with may others building
connections between the macro economy and asset value take aggregate
consumption as pre-specified. As we open the door to a richer collection
of macroeconomic models, it becomes important to entertain more
endogeneity, including investment and other variables familiar to
macroeconomics. In this section, we briefly describe one way to extend
the approach that builds directly on previous second-order approaches of
[@KimKimSchaumbergSims:2008], [@schmittgrohe_uribe:2004], and
[@LombardoUhlig:2018]. While such methods should not be viewed as being
generically applicable to nonlinear stochastic equilibrium models, we
find them useful pedagogically and often as at least initial steps to
understanding models that are arguably "smooth." See
[@PohlSchmeddersWilms:2018] for a careful study of nonlinearity in asset
pricing models with recursive utility.[^6]

We implement these methods for second-order approximation using the
following steps.

i)  Solve for ${\sf q}=0$ deterministic model.

ii) []{#stepi label="stepi"} Take as given a
    $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in $$\begin{aligned}
     \label{input_text} 
     {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 &=  {\frac 1 {1-\gamma_o}} \left[\mu^0 \cdot ( W_{t+1} - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \cr 
     {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2 & = {\frac 1 2} \left(W_{t+1} - \mu^0 \right)'\Upsilon_2^2 \left(W_{t+1} - \mu^0 \right) - 
     {\frac 1 2} \rm{tr}\left( \Upsilon_2^2\right) \cr & \hspace{.4cm} 
     +  \left(W_{t+1} - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right),
    \end{aligned}$$

iii) Compute the first-order expansion solve the resulting equations
     following the previous literature . When constructing these
     equations, use expectations computed using the probabilities
     induced by $N_{t+1}^0.$ Under the change in probability $W_{t+1}$
     is normally distributed with $\mu^0,$ where $\mu^0$ is given in
     step [\[stepi\]](#stepi){reference-type="ref" reference="stepi"}).
     Make an additional recursive utility adjustment to the equations
     that also depends on $\mu^0.$ It comes from a first-order
     approximation for $(\rho - 1)\left({\widehat V}_{t+1} 
     - {\widehat R}_t \right)$ as a plug in for the construction of the
     logarithm of the stochastic discount factor.

iv) Compute the second-order expansion and solve the resulting equations
    following the previous literature. Again use the expectations
    induced by $N_{t+1}^0$. In addition, make another recursive utility
    adjustment expressed in terms of $\mu^0, \Upsilon_0^2, \Upsilon_1^2$
    taking the inputs from [\[stepi\]](#stepi){reference-type="ref"
    reference="stepi"}. This comes from It comes from a second-order
    approximation for
    $(\rho - 1)\left({\widehat V}_{t+1} - {\widehat R}_t \right)$ needed
    for the logarithm of the stochastic discount factor.

v)  Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in
    representations [\[input_text\]](#input_text){reference-type="eqref"
    reference="input_text"} and return to
    [\[stepi\]](#stepi){reference-type="ref" reference="stepi"}). Repeat
    until convergence.

See Appendix [7](#appen:second_solve){reference-type="ref"
reference="appen:second_solve"} for more details and substantive
elaboration.

This first algorithm includes the use of first and second-order
approximation of $N_{t+1}^*.$ While the second-order approximation of
$N_{t+1}^*$ has conditional expectation one, it is not restricted to be
positive. A second algorithm enforces positivity on the implied
approximation of $N_{t+1}^*$.

i)  Solve for ${\sf q}=0$ deterministic model.

ii) []{#stepi2 label="stepi2"} Take as given an an exponential
    linear-quadratic approximation ${\widetilde N}_{t+1}$ for
    $N^*_{t+1},$ along with $\mu^0, \Upsilon_0^2, \Upsilon_1^2,$
    $\Upsilon_2^2$ used in formulas
    [\[input_text\]](#input_text){reference-type="eqref"
    reference="input_text"}.

iii) Compute the first-order expansion solve the resulting equations
     following the previous literature. When constructing these
     equations, use expectations computed using the probabilities
     induced by $N_{t+1}^a$. In addition, make an additional recursive
     utility adjustment to the equations that comes from a first-order
     approximation for $(\rho - 1)\left({\widehat V}_{t+1} 
     - {\widehat R}_t \right)$ as a plug in for the construction of the
     logarithm of the stochastic discount factor. Use the first-order
     adjustment from step [\[stepi2\]](#stepi2){reference-type="ref"
     reference="stepi2"}.

iv) Compute the second-order expansion and solve the resulting equations
    following the previous literature. Again use the expectations
    induced by ${\widetilde N}_{t+1}$. In addition, make another
    recursive utility adjustment expressed in terms of
    $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ taking the inputs from
    [\[stepi2\]](#stepi2){reference-type="ref" reference="stepi2"}. This
    comes from a second-order approximation for
    $(\rho - 1)\left({\widehat V}_{t+1}- {\widehat R}_t \right)$ needed
    for the construction of the logarithm of the stochastic discount
    factor.

v)  Form new values for
    $\mu^0, \Upsilon_0^2, \Upsilon_1^2, \Upsilon_2^2$ used in
    representations [\[input_text\]](#input_text){reference-type="eqref"
    reference="input_text"}, construct $${\widetilde  N}_{t+1}  =  
    \frac 
    {\exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right]} 
    {{\mathbb E}\left( \exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right] \mid {\mathfrak A}_t \right)},$$
    and return to [\[stepi2\]](#stepi2){reference-type="ref"
    reference="stepi2"}). Repeat until convergence.

vi) Compute

    $${ N}_{t+1}  =  
    \frac 
    {\exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right]} 
    {{\mathbb E}\left( \exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right] \mid {\mathfrak A}_t \right)} ,$$
    and set $N_{t+1}^a  = N_{t+1}$.

::: appendix
# Approximating solutions {#appen:second_solve}

For the purposes of this appendix, write:
$$\frac {S_{t+1}}{S_t} = N_{t+1}^* Q_{t+1} \beta \exp \left[- \rho \left( \log C_{t+1} - \log C_t \right) \right]$$
where $$\begin{aligned}
N_{t+1}^* & = \exp\left[  (1-\gamma_o) \left({\widetilde V}_{t+1} - {\widetilde  R}_t\right) \right] \cr
Q_{t+1}^* & = \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right]
\end{aligned}$$ are terms that are contributed by recursive utility.

## $N_{t+1}^*$ derivatives

$$\begin{aligned}
N_{t+1}^0 & \eqdef \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1}^0 - {\widetilde R}_t^0 \right) \right]   = \exp\left[  (1-\gamma_o) \left({\widehat V}_{t+1}^1 - {\widehat  R}_t^1\right) \right]\cr
N_{t+1}^1 & \eqdef \left.  \frac d {d{\sf q}}  \exp\left[(1 - \gamma_o) \left({\widetilde V}_{t+1} - {\widetilde R}_t \right) \right] \right\vert_{{\sf q} =0}  \cr & \hspace{.6cm} =  (1 - \gamma_o)  N_{t+1}^0\left( 
{\widetilde V}_{t+1}^1  - {\widetilde R}_t^1 \right) =   \left(\frac {1-\gamma_o} 2 \right) N_{t+1}^0 \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right).
\end{aligned}$$ It may be directly verified that $N_{t+1}^1$ has
conditional expectation equal to zero.

## $Q_{t+1}^*$ derivatives

$$\begin{aligned}
Q_{t+1}^0 & \eqdef \exp\left[(\rho - 1) \left({\widehat V}_{t+1}^0 - {\widehat R}_t^0 \right) \right]   = 1 \cr
Q_{t+1}^1 & \eqdef \left.  \frac d {d{\sf q}}  \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \right\vert_{{\sf q} =0}  = (\rho - 1) \left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \cr
Q_{t+1}^2  & \eqdef \left.  \frac {d^2}  {d{\sf q}^2}   \exp\left[(\rho - 1) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \right\vert_{{\sf q} =0}   = 
(\rho - 1)^2 \left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2 +  (\rho - 1) \left( 
{\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right)
\end{aligned}$$

Express $$\begin{aligned}
 \label{input_computation} 
 {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1 &=  {\frac 1 {1-\gamma_o}} \left[\mu^0 \cdot ( W_{t+1} - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \cr 
 {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2 & = {\frac 1 2} \left(W_{t+1} - \mu^0 \right)'\Upsilon_2^2 \left(W_{t+1} - \mu^0 \right) - 
 {\frac 1 2} \rm{tr}\left( \Upsilon_2^2\right) \cr & \hspace{.4cm} 
 +  \left(W_{t+1} - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right),
\end{aligned}$$ Recall from
[\[recur_update\]](#recur_update){reference-type="eqref"
reference="recur_update"} that
${\widehat V}_{t+1}^2 - {\widehat R}_{t}^2$ has mean zero under the
probability distribution induced by $N_{t+1}^0,$ which is consistent
with its representation in
[\[input_computation\]](#input_computation){reference-type="eqref"
reference="input_computation"}.

# Approximating expectation equations

Consider the equation:
$${\mathbb E} \left( N_{t+1} Q_{t+1} H_{t+1} \mid {\mathfrak A}_t \right)  + L_{t}   = 0 .$$
where
$\beta \exp \left[- \rho \left( \log C_{t+1} - \log C_t \right) \right]$
is absorbed into the construction of $H_{t+1}.$ This is the subsystem of
the equations not including the state evolution equations.

## Order zero

The order zero approximation of the product:
$N_{t+1} Q_{t+1} H_{t+1} +  L_{t}$ is:
$$N_{t+1}^0  H_{t+1}^0 + L_t^0 = 0$$ where we have substituted
$Q_{t+1}^0 = 1$. Thus the order zero approximate equation is:
$${\mathbb E} \left[N_{t+1}^0 \left( H_{t+1}^0 + L_{t+1}^0 \right)  \mid {\mathfrak A}_t \right] =  H_{t+1}^0 + L_{t}^0 = 0$$
since $N_{t+1}^0$ has conditional expectation equal to one. We add to
this subsystem the ${\sf q} = 0$ state dynamic equation inclusive of
jump variables, and we compute a stable steady state solution.

## Order one

The order one approximation of the product:
$N_{t+1} Q_{t+1} H_{t+1}+  L_{t}$ is:
$$N_{t+1}^1   H_{t+1}^0  + N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0 + N_{t+1}^0  H_{t+1}^1 +  L_{t}^1$$
where we have substituted $Q_{t+1}^0 = 1$. Thus the order one
approximate equation is: $$\begin{aligned}
& {\mathbb E} \left( N_{t+1}^1   H_{t+1}^0    + N_{t+1}^0 H_{t+1}^1  \mid {\mathfrak A}_t \right)
+ {\mathbb E} \left( N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right)  + L_t^1 \cr
& = {\mathbb E} \left[  N_{t+1}^0 \left(H_{t+1}^1  +  Q_{t+1}^1 H_{t+1}^0 \right) \mid {\mathfrak A}_t \right] + L_t^1 \cr
& = 0
\end{aligned}$$ where we used the implication that
$H_{t+1}^0  + L_{t+1}^0 = 0$. The contribution:
$${\mathbb E} \left(N_{t+1}^0  H_{t+1}^1   \mid {\mathfrak A}_t \right) + L_t^1$$
is of the form used for the first-order approximation without the
recursive utility modification, except that the expectation is evaluated
under the probability measure implied by $N_{t+1}^0$. The recursive
utility adjustment has us include the additional term:
$${\mathbb E} \left( N_{t+1}^0  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right) 
= {\frac {(\rho-1)}{2(1-\gamma_o)}} |\mu^o|^2H_{t+1}^0 \eqdef {\overline H}^1$$
which is constant over time. Thus we write the first-order subsystem of
equations as:
$${\mathbb E} \left(N_{t+1}^0   H_{t+1}^1 \mid {\mathfrak A}_t \right) + L_{t}^1 + {\overline H}^1 
 = 0.$$ We add to this the first-order approximation of the state
dynamics inclusive of jump variables and evaluate expectations under the
$N_{t+1}^0$ change of probability measure. Thus the one-period
conditional expectation of $W_{t+1}$ is $\mu^0.$

## Order two {#order-two}

The order two approximation of the product:
$N_{t+1} Q_{t+1} H_{t+1} + N_{t+1} L_{t+1}$ is: $$\begin{aligned}
& N_{t+1}^0 H_{t+1}^2 + N_{t+1}^2H_{t+1}^0  + 
2 N_{t+1}^1 H_{t+1}^1 + L_t^2 \cr + 
& 2 N_{t+1}^1  Q_{t+1}^1 H_{t+1}^0  +   2 N_{t+1}^0 Q_{t+1}^1 H_{t+1}^1 + N_{t+1}^0 Q_{t+1}^2 H_{t+1}^0 
\end{aligned}$$ The term $N_{t+1}^2H_{t+1}^0$ is zero and the term
$${\mathbb E} \left( N_{t+1}^0 H_{t+1}^2  \mid {\mathfrak A}_t \right) + L_t^2$$
coincides with the contribution for the second-order approximation
abstracting from recursive utility but evaluated under the change of
measure induced by $N_{t+1}^0$. Express $H_{t+1}^1$ as $$\begin{aligned}
 \label{first_solution}
H_{t+1}^1 & = \Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right). 
%L_{t+1}^1 & = \Lambda_0^1 +  \Lambda_1^1 X_t^1 + \Lambda_2^1 \left( W_{t+1} - \mu^0 \right) 
\end{aligned}$$ We now consider the additional terms $$\begin{aligned}
 \label{second_affine} 
2{\mathbb E}\left( N_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 
(1 - \gamma_o)  {\mathbb E} \left[ N_{t+1}^0   \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right)H_{t+1}^1 \mid{\mathfrak A}_t \right],\cr 
& = (1 - \gamma_o)  \Theta_2^1 \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2 \right) \cr \cr
 2 {\mathbb E} \left(N_{t+1}^1  Q_{t+1}^1 H_{t+1}^0 \mid {\mathfrak A}_t \right) & = (\rho  - 1)(1 - \gamma_o){\mathbb E}\left[ N_{t+1}^0 \left( {\widehat V}_{t+1}^2 - {\widehat R}_{t}^2\right) \left( {\widehat V}_{t+1}^1 - {\widehat R}_{t}^1\right)  \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
 & = (\rho  - 1) \mu^o \cdot \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2 \right) H_{t+1}^0\cr \cr
2 {\mathbb E} \left(N_{t+1}^0 Q_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 2(\rho - 1) {\mathbb E}\left[ N_{t+1}^0\left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right) \right] \mid {\mathfrak A}_t \right] \cr
&  = 2 \frac {(\rho - 1)}{(1-\gamma_o)}\left[   \Theta_2^1 \mu^0 + {\frac 1 2} \mu^0 \cdot \mu^0 \left(\Theta_0^1 + \Theta_1^1 X_t^1\right) \right] \cr \cr
{\mathbb E} \left(  N_{t+1}^0 Q_{t+1}^2 H_{t+1}^0 \mid {\mathfrak A}_t \right) & =  (\rho - 1)^2{\mathbb E} \left[N_{t+1}^0 \left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2  \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr & = \left(\frac {1-\rho}{1-\gamma_o}\right)^2\left(|\mu^0|^2 + {\frac 1 4} |\mu^0|^4 \right) H_{t+1}^0
\end{aligned}$$ Denote the sum of the four terms in
[\[second_affine\]](#second_affine){reference-type="eqref"
reference="second_affine"} as ${\overline H}_t^2.$ This random variable
will be affine in $X_t^1$, with a dynamic evolution determined by
solving the first-order approximation. Thus we write the subsystem of
equations to be solved as:
$${\mathbb E} \left( N_{t+1}^0 H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\overline H}_t^2  = 0.$$
We add to this second-order subsystem, the second-order approximation of
the state dynamics inclusive of the of jump variables. We substitute in
the solution for the first-order approximation for the jump variables
into both the first and second-order approximate state dynamics. In
solving the second-order jump variable adjustment we use expectations
induced by $N_{t+1}^0$ zero throughout under which $W_{t+1}$ is
conditionally normally distributed with mean $\mu^0$ and covariance $I$.

# Steps for implementation

We implement these methods for second-order approximation using the
following steps.

i)  Solve $H_{t+1}^0$ and $L_{t+1}^0$ for order zero state and jump
    variables. The outcome will be state invariant.

ii) []{#stepi label="stepi"} Take as given a
    $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in representations
    [\[input_computation\]](#input_computation){reference-type="eqref"
    reference="input_computation"}.

iii) Compute the first-order contribution to approximation by following
     the previous literature with expectations computed using the
     probabilities induced by $N_{t+}^0$, which imply that $W_{t+1}$ has
     mean $\mu^0.$ Express the solution as in
     [\[first_solution\]](#first_solution){reference-type="eqref"
     reference="first_solution"}.

iv) Compute the second-order contribution to the approximation by
    following the previous literature, again with the expectations
    induced by $N_{t+1}^0.$

v)  Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2$ used in
    representations
    [\[input_computation\]](#input_computation){reference-type="eqref"
    reference="input_computation"} and return to
    [\[stepi\]](#stepi){reference-type="ref" reference="stepi"}). Repeat
    until convergence.

# Second approach

For this solution, we iterate over $N_{t+1}^*$ approximation. Call the
approximation ${\widetilde N}_{t+1}$ with an induced distribution for
$W_{t+1}$ that is normal with conditional mean ${\tilde \mu}_t$ and
covariance matrix ${\widetilde \Sigma}$. This distribution is used in
both the first-order and second-order contributions to the approximation
The conditional mean for ${\tilde \mu}_t$ is affine in $X_t^1$. The
following delineates the changes that need to be made.

## First-order adjustment

Compute: $$\begin{aligned}
{\mathbb E} \left( {\widetilde N}_{t+1}  Q_{t+1}^1 H_{t+1}^0  \mid {\mathfrak A}_t \right) & = 
(\rho - 1) 
{\mathbb E} \left[ {\widetilde N}_{t+1} \left({\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
& = \left({\frac {\rho -1} {1-\gamma_o}}\right)  \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] H_{t+1}^0 \cr 
& \eqdef {\widetilde  H}_t^1.
\end{aligned}$$ Then the equation to be solved is:
$${\mathbb E} \left({\widetilde N}_{t+1}   H_{t+1}^1 \mid {\mathfrak A}_t \right) + L_{t}^1 + {\widetilde H}^1_t  
 = 0.$$

## Second-order adjustment

$$\begin{aligned}
 \label{second_affine_again} 
2 {\mathbb E} \left({\widetilde N}_{t+1} Q_{t+1}^1 H_{t+1}^1 \mid {\mathfrak A}_t \right) & = 2(\rho - 1) {\mathbb E}\left[ {\widetilde N}_{t+1}\left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right) \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1 \left( W_{t+1} - \mu^0 \right) \right] \mid {\mathfrak A}_t \right] \cr
&  = 2 \frac {(\rho - 1)}{(1-\gamma_o)}  \Theta_2^1 {\widetilde \Sigma} \mu^0 \cr 
& \hspace{.5cm} + 2 \frac {(\rho - 1)}{(1-\gamma_o)} \left[\mu^0 \cdot ( {\tilde \mu}_t - \mu^0) +  {\frac 1 2} \mu^0 \cdot \mu^0\right] \left[\Theta_0^1 + \Theta_1^1 X_t^1 + \Theta_2^1\left( {\tilde \mu}_t - \mu^0\right)  \right] 
\cr \cr
{\mathbb E} \left(  {\widetilde N}_{t+1} Q_{t+1}^2 H_{t+1}^0 \mid {\mathfrak A}_t \right) & =  (\rho - 1)^2{\mathbb E} \left[{\widetilde N}_{t+1} \left( 
{\widehat V}_{t+1}^1  - {\widehat R}_t^1 \right)^2  \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr 
& \hspace{.5cm} +   (\rho - 1) {\mathbb E} \left[{\widetilde N}_{t+1}  \left( 
{\widehat V}_{t+1}^2  - {\widehat R}_t^2 \right) \mid {\mathfrak A}_t \right] H_{t+1}^0 \cr
& = \left(\frac {1-\rho}{1-\gamma_o}\right)^2\left[{\mu^0}'{\widetilde \Sigma} \mu^0 + \left(\mu^0 \cdot {\tilde \mu}_t - {\frac 1 2} |\mu^0|^2 \right)^2\right]
H_{t+1}^0 \cr
& \hspace{.5cm} + \frac {(\rho - 1)} 2 \left[ \rm{tr}\left( \Upsilon_2^2{\widetilde \Sigma} - \Upsilon_2^2 \right) +  \left({\tilde \mu}_t - \mu^0\right)' \Upsilon_{2}^2  \left({\tilde \mu}_t - \mu^0\right)\right] H_{t+1}^0 \cr
& \hspace{.5cm} +  (\rho - 1) \left({\tilde \mu}_t  - \mu^0\right)' \left( \Upsilon_1^2 X_t^1 + \Upsilon_0^2\right) H_{t+1}^0
\end{aligned}$$ Denote the sum of the two terms in
[\[second_affine_again\]](#second_affine_again){reference-type="eqref"
reference="second_affine_again"} as ${\widetilde  H}_t^2.$ Then the
equation to be solved is
$${\mathbb E} \left( {\widetilde N}_{t+1} H_{t+1}^2 \mid {\mathfrak A}_t \right) + L_{t}^2 + {\widetilde H}_t^2  = 0.$$

## Updated recursive utility adjustments

Form new values for $\mu^0, \Upsilon_0^2, \Upsilon_1^2,  \Upsilon_2^2$
used in representations
[\[input_computation\]](#input_computation){reference-type="eqref"
reference="input_computation"}. Compute a new version of
$${\widetilde  N}_{t+1}  =  
\frac 
{\exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right]} 
{{\mathbb E}\left( \exp\left[ (1 - \gamma_o) \left[ {\widehat V}^1_{t+1} -  {\widehat R}^1_{t}+ \frac 1  2 \left({\widehat V}^2_{t+1} -{\widehat R}^2_{t}  \right) \right] \right] \mid {\mathfrak A}_t \right)},$$
and deduce the implied ${\widetilde \mu}_t$ and ${\widetilde \Sigma}$.
The conditional mean ${\tilde \mu}_t$ satisfies:
$${\widetilde \Sigma}^{-1} {\tilde \mu}_t =  \mu^0 +  {\frac {(1-\gamma_o)} 2} \left(\Upsilon_0^2 + \Upsilon_1^2 X_t^1 - \Upsilon_2^2 \mu^0 \right)$$
where the formula for ${\widetilde \Sigma}$ is
$$\widetilde {\Sigma} =  \left[{\mathbb I} - {\frac {(1 - \gamma_o)}  2} \Upsilon_2^2\right]^{-1}.$$

With these adjustments, we iterate to convergence.

# Parameter values

    parameter                              value                                       source
  ------------- ------------------------------------------------------------ ---------------------------
   $\sigma_k$     $.01 \times \begin{bmatrix} .477 & 0 & 0 \end{bmatrix}$        HS tenuous beliefs
   $\sigma_z$    $.01 \times \begin{bmatrix} .011 & .025 & 0 \end{bmatrix}$      HS tenuous beliefs
   ${\sf a}_z$                              .986                                 HS tenuous beliefs
   $\sigma_y$          $\begin{bmatrix} 0 & 0 & .0648 \end{bmatrix}$          SSY Econometrica, table 3
   ${\sf a}_y$                              .964                              SSY Econometrica, table 7
   $\sigma_y$           $\begin{bmatrix} 0 & 0 & .108 \end{bmatrix}$          SSY Econometrica, table 3
   ${\sf a}_y$                              .976                              SSY Econometrica, table 7
    $\delta$                                .002                                    HBSH twisted

  : Parameters used for quarterly models.

The direct shock to capital is ${\sf q} \exp(Y_t) \sigma_k W_{t+1}$.
$$\begin{aligned}
Z_{t+1}  =  & {\sf a}_z Z_t + {\sf q} \exp(Y_t) \sigma_z W_{t+1} \cr
Y_{t+1} = & {\sf a}_y Y_t + {\sf q} \sigma_y W_{t+1} \cr
\log K_{t+1} - \log K_t = & \left(\frac {I_t}{K_t} \right) - 14  \left(\frac {I_t}{K_t} \right)^2 
- .0128 + Z_t  - \frac {{\sf q}^2  \exp(2Y_t) |\sigma_k|^2} 2  + {\sf q} \exp(Y_t) \sigma_k W_{t+1}
\end{aligned}$$ Notice that the Jensen term only contributes a constant
to the second order.

Reconfigure the shock process as follows. Form
$$\begin{bmatrix} .477 & 0 \cr .011 & .025 \end{bmatrix} \begin{bmatrix} .477 & .011 \cr 0  & .025 \end{bmatrix}  =
\begin{bmatrix} ? & ? \cr 0 & ? \end{bmatrix} \begin {bmatrix} ?  & 0 \cr ?  & ?  \end{bmatrix}$$
Take $\sigma_k$ to be the first row of left side and $\sigma_z$ the
second row of the left side.

For the endowment economy use
$$\log C_{t+1} - \log C_t = .00484 + Z_t + \exp(Y_t) \sigma_k W_{t+1}$$
:::

[^1]: Some background notes for the Macro International Workshop at the
    University of Chicago. The approximations and computations described
    in these notes are supported by a jupyter notebook entitled
    'uncertainexpansion.jpynb' in a repository
    https://github.com/lphansen/RiskUncertaintyValue. This repository
    also contains a complementary jupyter notebook entitled
    'shockelasticity.jpynb' that computes impulse responses and shock
    elasticities for models represented by their second-order
    approximations.

[^2]: [@LombardoUhlig:2018] provide a discussion of how their approach
    builds on more general perturbation methods as discussed by
    [@Holmes:2012] and [@Judd:1998].

[^3]: The [@HansenHeatonLi:2008] predictability evidence turned out to
    be "fragile" and was modified and updated in [@HansenSargent:2021]
    Appendix B. This same appendix suggests a way to deduce a
    statistical approximation to the first order dynamics of
    [@bansalyaron2004] from a more general VAR representation of the
    consumption dynamics.

[^4]: It is notable that we are looking at levels and not logarithms of
    consumption. the local impulse response for the logarithms of
    consumption is in fact zero for the stochastic volatility shock.

[^5]: We normalized the stochastic volatility shock $\sigma_x^2$ to be
    negative implying that a positive shock reduces the stochastic
    volatility state variable. Under this normalization, the shock price
    elasticities are positive.

[^6]: @PohlSchmeddersWilms:2018 provide examples of when log-linear or
    local methods of computation fail to provide good approximations.
