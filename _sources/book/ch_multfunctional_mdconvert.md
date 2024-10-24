(chap:mult)=
# Multiplicative Functionals

Chapter [](chap:add) described additive functionals of a Markov process. This chapter describes exponentials of additive functionals that we call *multiplicative functionals* and that we can use to model stochastic growth, stochastic discounting, and their interactions. After adjusting for geometric growth or decay, a multiplicative functional contains a martingale component that turns out to be a likelihood ratio process that is itself a special type of multiplicative functional called an exponential martingale. In chapter [](chap:like), we shall see how a likelihood ratio process is used to construct an alternative probability model with Markov dynamics by simply multiplying an original probability measure by the likelihood ratio process. This procedure is widely used in asset pricing theory and helps us to represent stochastic components of growth and discounting that persist over long horizons. To analyze multiplicative functionals, we apply mathematical tools closely related to tools from the statistical theory of large deviations that we shall apply again in chapter [](chap:like).

## Geometric Growth and Decay

To construct a multiplicative functional, we start with an underlying Markov process $X$ that has stationary distribution $Q$.

````{prf:definition}
:label: def:multfunctional
Let $\{ Y_t\}$ be an additive functional that as in Chapter [](chap:add) is described by $Y_{t+1} - Y_t = \kappa(X_t, W_{t+1})$, where $X_t$ is the time $t$ component of a Markov state vector and $W_{t+1}$ is the time $t+1$ value of an i.i.d. process of unanticipated shocks. We say that $\{M_t \} = \{ \exp(Y_t) \}$ is a **multiplicative functional** parameterized by $\kappa$. When $Y_0$ is a (Borel measurable) function of $X_0$, $M_0 >0$ is also a (Borel measurable) function of $X_0$.
````

An additive functional grows or decays *linearly*, so the exponential of an additive functional grows or decays *geometrically*. In Chapter [](chap:add) we constructed a Law of Large Numbers and a Central Limit Theorem for additive functionals. In this chapter, we use other mathematical tools to analyze the limiting behavior of multiplicative functionals.

(sec:threespecial)= 
## Special Multiplicative Functionals

We define the three primitive types of multiplicative functionals.


````{prf:example}
:label: exmult1
Suppose that $\kappa = \bar \kappa$ is constant and that $M_0$ is a Borel measurable function of $X_0$. Then
```{math}
M_t  = \exp\left( t \bar \kappa \right)M_0.
```
This process grows or decays geometrically.
````


````{prf:example}
:label: exmult2
Suppose that
```{math}
E\left[ \exp \left[ \kappa \left( X_t, W_{t+1} \right) \right] \vert X_t \right] = 1.
```
Then
```{math}
:label: eqn:multcomp2
E\left(M_{t+1} \vert \mathcal{F}_t \right) = M_t .
```
A multiplicative functional that satisfies {eq}`eqn:multcomp2` is called a *multiplicative martingale*.
````


````{prf:example}
:label: exmult3
Suppose that $M_t = \exp\left[h(X_t)\right]$ where $h$ is a Borel measurable function. The associated additive functional satisfies
```{math}
Y_{t+1} - Y_t = \log M_{t+1} - \log M_t = h(X_{t+1}) - h(X_t) = h\left[ \phi(X_t, W_{t+1} ) \right] - h(X_t)
```
and so is parameterized by $\kappa(X_t, W_{t+1}) = h\left[ \phi(X_t, W_{t+1} ) \right] - h(X_t)$ with initial condition $Y_0 = h(X_0)$.
````

When the process $\{X_t\}$ is stationary and ergodic, multiplicative functional {prf:ref}`exmult1` displays expected growth or decay, while multiplicative functionals {prf:ref}`exmult2` and {prf:ref}`exmult3` do not. Multiplicative functional {prf:ref}`exmult3` is stationary, while {prf:ref}`exmult1` and {prf:ref}`exmult2` are not.

We can construct other multiplicative functionals simply by multiplying two or more instances of these primitive ones. In the following section, we reverse that process: we take an arbitrary multiplicative functional and (multiplicatively) decompose it into instances of our three types of multiplicative functionals.


## Multiplicative martingales and likelihood processes

Multiplicative martingales induce alternative probabilities and are thus likelihood-ratio processes.  



(sec:factorization)=
## Factoring Multiplicative functionals

Following {cite}`hansenscheinkman09` and {cite}`hansen:2012`, we decompose a multiplicative functional into multiplicative components having the primitive types {prf:ref}`exmult1`, {prf:ref}`exmult2`, {prf:ref}`exmult3`. As in definition {prf:ref}`def:multfunctional`, let $\{Y_t\}$ be an additive functional and let $M_t = \exp(M_t)$. To bounded Borel measurable functions $f$ of the Markov state, we apply a one-period operator $\mathbb{M}$ defined by

```{math}
:label: Mdef
\mathbb{M}f(x) = E\left[ \exp(Y_{t+1} - Y_t) f(X_{t+1}) \vert X_t = x \right] = E\left[ \frac{M_{t+1}}{M_t} f(X_{t+1}) \vert X_t = x \right].
```

A two-period operator is

```{math}
:label: eq:M2def
\mathbb{M}^2 f(x) = E\left[ \exp(Y_{t+2} - Y_t) f(X_{t+2}) \vert X_t = x \right] = E\left[ \frac{M_{t+2}}{M_t} f(X_{t+2}) \vert X_t = x \right],
```

with corresponding definitions of $j$-period operators $\mathbb{M}^j$. For $f > 0$, take the limiting growth rate to be

```{math}
:label: eqn:limgrowthrate
\eta = \lim_{j\rightarrow \infty} \frac{1}{j} \log E\left[ \exp(Y_{t+j} - Y_t) f(X_{t+j}) \bigl | X_t = x \right] = \lim_{j\rightarrow \infty} \frac{1}{j} \log E\left[ \frac{M_{t+j}}{M_t} f(X_{t+j}) \bigl| X_t = x \right] = \lim_{j\rightarrow \infty} \frac{1}{j} \log \mathbb{M}^j f(x).
```

Taking $f$ to be unity implies that

```{math}
:label: eqn:etadef
\eta = \lim_{j\rightarrow \infty} \frac{1}{j} \log E\left[ \frac{M_{t+j}}{M_t} \Bigl| X_t = x \right].
```

We call $\eta$ the expected asymptotic growth rate of the multiplicative functional $\{M_t\}_{t=0}^\infty$. Multiplying the multiplicative functional by $\exp(-\eta t)$ removes expected asymptotic growth.

To refine our characterization of a multiplicative functional, we pose:

(prob:eigen)=
**Eigenvalue-eigenfunction Problem:**
Solve

```{math}
:label: eqn:eigenmult
\mathbb{M}{\tilde e}(x) = \exp\left( {\tilde \eta} \right) {\tilde e}(x)
```
for an eigenvalue $\exp(\tilde \eta)$ and a positive eigenfunction ${\tilde e}$.


We call the largest eigenvalue the principal eigenvalue and the associated eigenvector the principal eigenfunction of the operator $\mathbb{M}$.

A positive eigenfunction ${\tilde e}$ is a function of the Markov state that can be expected to grow geometrically at the long-run growth rate $\eta = {\tilde \eta}$. Write the eigenfunction equation {eq}`eqn:eigenmult` as:

[^teachingNote]: Iterating the eigenfunction equation implies $E\left[ {\frac {{ M}_{t+j}}{{ M}_t}} {\tilde e}(X_{t+j}) \vert X_t \right] =  \exp\left( { j \tilde \eta} \right) {\tilde e}(X_t).$ 

```{math}
E\left[ {\frac {{ M}_{t+1}}{{ M}_t}} {\tilde e}(X_{t+1}) \vert X_t \right] =  \exp\left( {\tilde \eta} \right) {\tilde e}(X_t).
```

Solve for the principal eigenvalue and eigenvector. Then define

```{math}
{\frac {{\widetilde M}_{t+1}}{{\widetilde M}_t}} \equiv \exp(- {\tilde \eta}) {\frac {M_{t+1}{\tilde e}(X_{t+1})}{M_t {\tilde e}(X_t)}} = \exp\left[{\tilde \kappa}(X_t, W_{t+1}) \right].
```

Evidently, ${\frac {{\widetilde M}_{t+1}}{{\widetilde M}_t}}$ has a conditional expectation equal to unity. Consequently, $\{ {\widetilde M}_t  \}$ is a multiplicative martingale with $\widetilde M_0 = 1$.

````{prf:theorem}
:label: prop:factor

Let $\{ M_t \}$ be a multiplicative functional.
Suppose that the principal eigenvalue-eigenfunction [Problem](prob:eigen) has a solution with principal eigenfunction
${\tilde e}(X)$.
Then the multiplicative functional
is the product of three components that are instances of the primitive functionals in examples {prf:ref}`exmult1`, {prf:ref}`exmult2`, and {prf:ref}`exmult3`:
```{math}
:label: eqn:multdecomposition
{\frac {M_t}{ M_0}} = \exp\left( {\tilde \eta} t \right) \left[ {\frac {{\widetilde M}_t}{{\widetilde M}_0 }} \right] \left[ {\frac {{\tilde e}(X_0)}{{\tilde e}(X_t) }} \right]
```
where $\{ {\widetilde M}_t \}$ is a multiplicative martingale.
````

The decomposition of a multiplicative functional described in {prf:ref}`prop:factor` is a counterpart to the {prf:ref}`prop:decomp` decomposition of an additive functional.

````{prf:remark}
If there is a martingale component in an additive (in logarithms) decomposition, then there is necessarily a martingale in the corresponding multiplicative (in levels) {prf:ref}`prop:factor` factorization, and conversely.
````

````{prf:theorem}
Because $\widetilde M_0 = 1$, the multiplicative martingale $\widetilde M_t$ is a likelihood ratio process.
````

````{prf:remark}
We used the martingale in {prf:ref}`prop:decomp` to identify the permanent component of an additive functional in Chapter [](chap:add). In this chapter, we shall use the multiplicative martingale isolated by {prf:ref}`prop:factor` to represent a change of probability measure.
````

The fact that the additive martingale $Y_t = \log(M_t)$ has a variance that grows linearly over time contributes a component to the exponential trend of the multiplicative functional $\{M_t\}$. The following log-linear, log-normal model displays relevant mechanics.

````{prf:example}
:label: ex:VARmultfunctional

Consider a stationary $X$ process described by the VAR

```{math}
X_{t+1} = A X_t + B W_{t+1},
```
where $A$ is a stable matrix and
$\{ W_{t+1} \}$ is a sequence of independent and identically
normally distributed random vectors with mean zero and covariance matrix $I$.
Let $\{Y_t\}$ be the additive functional generated by

```{math}
Y_{t+1} - Y_t = \kappa(X_{t},W_{t+1}) = \nu + D \cdot X_t + F \cdot W_{t+1},
```
where $D$ and $F$ are vectors with the same dimensions as $X_t$ and $W_{t+1}$, respectively.
In Proposition {prf:ref}`prop:decomp` of Chapter [](chap:add), we described the decomposition

```{math}
:label: adddecomp2
Y_{t} - Y_0 =  t\nu + \left[\sum_{j=1}^{t} \kappa_a(X_{j-1},W_{j})\right]   - g(X_{t}) + g(X_0)
```
where
$\kappa_a(X_{t},W_{t+1}) = H \cdot W_{t+1}, H =  [F + B'(I-A')^{-1} D]$,
and
$g(x) = D'(I-A)^{-1} x$.
Let $M_t = \exp(Y_t)$. Use equation {eq}`adddecomp2` to deduce

```{math}
\frac{M_t}{M_0} = \exp (t \nu) \exp \Bigl(\sum_{j=1}^t H \cdot W_j \Bigr) \exp \biggl( D'(I-A)^{-1} X_0 - D'(I-A)^{-1} X_t \biggr)
```
or

```{math}
\frac{M_t}{M_0} =  \exp\left( \tilde \eta t \right) \Biggl( \frac{\widetilde M_t}{\widetilde M_0}\Biggr) \left( \frac{\tilde e (X_0)} {\tilde e(X_t)} \right)
```
where

```{math}
\tilde \eta =  \nu + \frac{H \cdot H}{2} ,
```

```{math}
\widetilde M_t = \exp \biggl( \sum_{j=1}^t \biggl(H \cdot W_j -\frac{ H \cdot H }{2} \biggr) \biggr),  \quad \widetilde M_0 =1 ,
```

and

```{math}
\tilde e(x) = \exp[g(x)] = \exp \bigl[ D' (I - A)^{-1} x \bigr].
```

````

````{prf:example}
:label: ex:tallarini

{cite}`tallarini2000` compared utility functionals associated with two models of per capita US consumption, each being a multiplicative process taking the form

```{math}
:label: eqn:tall101
c_{t+1} = \mu + c_t + \sigma_{\epsilon} \epsilon_{t+1},
```

but having different $(\mu, \sigma_\epsilon)$ parameter values, where $\{\epsilon_{t+1}\}_{t=0}^\infty$ is an i.i.d. sequence of $\mathcal{N}(0,1)$ random variables, $c_t = \log C_t$, and $C_t$ is per capita consumption of nondurables and services. {cite}`tallarini2000` chose this random walk plus drift specification for log per capita consumption because it approximates post World War II U.S. data well. Consumption process {eq}`eqn:tall101` is an additive functional that is a special case of our VAR specification with $A=0, B=0, D = 0, F = \sigma_{\epsilon}$. Process {eq}`eqn:tall101` for $c_t$ implies that

```{math}
C_{t+1} = C_t \exp(\mu) \exp(\sigma \epsilon_{t+1}) = C_t \exp\left(\mu + \frac{\sigma^2}{2} \right) m(\epsilon_{t+1})
```

where $m(\epsilon_{t+1}) = \exp\left( \sigma \epsilon_{t+1} - \frac{\sigma^2}{2} \right) \geq 0$ and $E m(\epsilon_{t+1} ) =1$, so that $m(\epsilon_{t+1})$ is a likelihood ratio and

```{math}
E C_{t+1} | C_t = \exp\left( \mu + \frac{\sigma^2}{2} \right) C_t.
```

The level of consumption is the stochastic process

```{math}
:label: eqn:tall102
C_t = C_0\exp\left(\mu + \frac{\sigma_{\epsilon}^2}{2}\right) \exp\left( \sum_{j=1}^t \left(\sigma_\epsilon \epsilon_j - \frac{\sigma_{\epsilon}^2}{2}\right) \right),
```

where $\exp\left( \sum_{j=1}^t \left(\sigma_\epsilon \epsilon_j - \frac{\sigma_{\epsilon}^2}{2}\right) \right) = \prod_{j=1}^t m(\epsilon_j)$ is a multiplicative martingale with mean 1.[^likelihoodratio]

{cite}`tallarini2000` compared the utility consequences of two processes (or "consumption plans") with identical exponential trends $\exp\left( \mu + \frac{\sigma^2}{2} \right)$ but differing shock processes $m(\epsilon_{t+1})$, one being a risky plan with $m(\epsilon_{t+1}) = \exp\left( \sigma \epsilon_{t+1} - \frac{\sigma^2}{2} \right)$ being a nontrivial random variable, the other being a no-risk process with $m(\epsilon_{t+1}) \equiv 1$ for all $t \geq 0$. {cite}`tallarini2000` calibrated the risky plan A to match quarterly post World War II U.S. per capita aggregate consumption by setting $\sigma_{\epsilon} = \sigma_{\epsilon}^A \approx .005$ and $\mu = \mu^A \approx .005$. Plan B is a hypothetical "no-risk" path that (i) sets $\sigma_\epsilon = \sigma_\epsilon^B =0$ and $\mu = \mu^B = \mu^A + \frac{(\sigma_{\epsilon}^{A})^2}{2}$, and (ii) adjusts $C_0$ downward to compensate for the reduced risk in plan B relative to plan A. According to equation {eq}`eqn:tall102`, the two plans have the same geometric growth rate $\mu^A + \frac{(\sigma_{\epsilon}^{A})^2}{2}$. But under plan $A$, the process $\{C_t\}_{t=0}^\infty$ also has risk in the form of the multiplicative martingale $\exp\left( \sum_{j=1}^t \left(\sigma_\epsilon^A \epsilon_j - \frac{(\sigma_{\epsilon}^{A})^2}{2}\right) \right)$, a component that is absent from the no-risk plan B. Tallarini assumed that a representative consumer evaluates these alternative consumption plans $\{C_t\}_{t=0}^\infty$ according to the following time $0$ utility criterion that is the $A = B = D =0, F = \sigma_\epsilon$ special case of the Proposition {prf:proposition}`prop:tallvalue` risk-sensitive utility recursion

```{math}
:label: eqn:tall103
U_0 = \frac{\beta}{1-\beta} \left[ \mu - \frac{\sigma_\epsilon^2}{2 \theta (1 - \beta)} \right] + c_0,
```

where $\theta = \frac{-1}{(1-\beta)(1 - \gamma)}$, $\beta \in (0,1)$ is a discount factor, $\gamma \geq 1$ is a coefficient of relative risk aversion, and $\sigma_\epsilon$ equals zero under plan B and $\sigma_\epsilon = \sigma^A \approx .005$ under plan A. {cite}`tallarini2000` calculated a percentage reduction in $c_0^B$, the logarithm of initial plan B consumption, that would make the representative consumer indifferent between the risky plan A and the no-risk plan B. Utility indexes $U_0$ under the two plans are

```{math}
U_0^A & = \frac{\beta}{1 - \beta} \left[ \mu_A - \frac{\sigma_\epsilon^2}{2 \theta (1 - \beta)} \right] + c_0^A \\
U_0^B & = \frac{\beta}{1 - \beta} \left[ \mu_A + \frac{\sigma_\epsilon^2}{2} \right] + c_0^B
```

Taking $c_0^A$ as a fixed value of $c_0$ under plan $A$, equating $U_0^A$ to $U_0^B$, and solving for $c_0^B $ gives

```{math}
:label: eqn:tallwelfcosts
c_0^B - c_0^A = \frac{\beta}{1-\beta} \left[ \frac{- (\sigma_\epsilon^A)^2 \gamma }{2} \right].
```

Formula {eq}`eqn:tallwelfcosts` tells us the percentage reduction in time $0$ consumption that a representative consumer would be willing to accept in order not to bear the risk associated with the multiplicative martingale component of plan A.

[^likelihoodratio]: Therefore, it is a likelihood ratio process.
````

````{prf:remark}
Formula {eq}`eqn:tall103` is a special case of the time zero valuation function described in Proposition {prf:proposition}`prop:tallvalue` that is appropriate for the special additive process for $\{\log C_t\}_{t=0}^\infty$ assumed by {cite}`tallarini2000`.
````

````{prf:theorem}
(Forecasting a Multiplicative Functional, I)

Let $\{Y_t\}$ be the additive functional
```{math}
:label: eqn:addex1b2
Y_{t+1} - Y_t = \kappa(X_{t},W_{t+1}) = \nu + D \cdot X_t + F \cdot W_{t+1},
```
where $D$ and $F$ are vectors with the same dimensions as $X_t$ and $W_{t+1}$, respectively. It follows that
```{math}
\log Y_{t+j} - \log Y_t = j \nu + D (I + A + \cdots A^{j-1}) X_t + F W_{t+j}
+ (F + DB) W_{t+j-1} + (F + D(I+A)B) W_{t+j-2}
+ \cdots + (F + D(I + A + \cdots + A^{j-2}) B ) W_{t+1}
```
which it is convenient to express as
```{math}
\log Y_{t+j} - \log Y_t = j \nu + \check G_j X_t + \sum_{k=0}^{t-1} \check H_k W_{t+j-k},
```
where $\check G_j$ and $\check H_j$ are defined to make this expression match the previous one. Where
```{math}
\Omega_j =\sum_{k=0}^{j-1} H_k \cdot H_k,
```
it follows that
```{math}
\log Y_{t+j} - \log Y_t \sim {\mathcal N}(j \nu + \check G_j X_t, \Omega_j)
```
and that
```{math}
:label: eqn:Emultfunct
E \left[ \frac{Y_{t+j}}{Y_t} \right]\Bigl| {\mathfrak F}_t = \exp\left( j v + \check G_j X_t + \frac{1}{2} \Omega_j\right).
```
````

````{prf:example}
(Term Structure of Interest Rates, I)

Let the Markov state $X_t$ be governed by the VAR
```{math}
X_{t+1} = A X_t + B W_{t+1},
```
where $A$ is a stable matrix and
$\{ W_{t+1} \}$ is a sequence of independent and identically
normally distributed random vectors with mean zero and covariance matrix $I$.
Let consumption $\{C_t\}_{t=0}^\infty$ be the additive functional generated by
```{math}
C_{t+1} - C_t =  \nu + D \cdot X_t + F \cdot W_{t+1}, \quad t \geq 0,
```
where $D$ and $F$ are vectors with the same dimensions as $X_t$ and $W_{t+1}$, respectively.
A consumer ranks consumption streams $\{C_t\}_{t=0}^\infty$
according to
```{math}
E \left[ \sum_{t=0}^\infty \exp(-\delta t) \frac{C_t^{1-\gamma}}{1-\gamma} \right] \Bigg| { \mathfrak F}_0 ,
```
where $\delta >0$ is a rate of time preference and $\gamma \geq 1$ is a coefficient of relative risk aversion.
This consumer has a stochastic discount factor process $\{S_t \}_{t=0}^\infty$ governed by
```{math}
\log S_{t+1} - \log S_t = -\delta - \gamma (\log C_{t+1} - \log C_t )
```
or
```{math}
\log S_{t+1} - \log S_t = - (\delta + \gamma \nu) - \gamma(D \cdot X_t + F \cdot W_{t+1}) .
```
We can apply formula {eq}`eqn:Emultfunct` to compute the price of a claim at time $t$ to a risk-free claim
on one unit of consumption at time $t+1$
```{math}
E \left[ \frac{S_{t+j}}{S_t} \right]\Big| {\mathfrak F}_t ,
```
thereby acquiring a theory of the term structure of interest rates.
````

### Finite State Markov Increments

The case in which $X_t$ is a finite-state Markov chain is also manageable computationally.

````{prf:example}
:label: ex:Markov_chain_mult
(Finite state Markov chain)

The stationary (or asymptotically stationary) stochastic process $X_t$ is governed by a finite state Markov chain on state space
${\mathbf S} = \{ s_1, s_2, \ldots, s_n \}$, where
$s_i$ is the $n \times 1$ unit vector whose components are all zero except for $1$ in the $i$th row.  The transition matrix
is $P$ where $P_{ij} = \textrm{Prob}( X_{t+1} = s_j | X_t = s_i)$.  We can represent the Markov chain as
```{math}
X_{t+1} = P' X_t + W_{t+1}
```
where $E (X_{t+1} | X_t ) = P' X_t $, $P'$ denotes the transpose of $P$,  and $\{W_{t+1}\}$ is an $n \times 1$ vector process that satisfies $E ( W_{t+1} | X_t) = 0 $, which is therefore a
martingale difference sequence adapted to $X_t, X_{t-1}, \ldots , X_0$.

Let $G$ be an $n \times n$ matrix whose $(i,j)$ entry $G_{ij}$ is an additive net  growth rate  $Y_{t+1} - Y_t$ experienced   when $X_{t+1} = s_j$ and
$X_t = s_i$.  The stochastic process $\{Y_t\}$ is governed by the *additive functional*
```{math}
Y_{t+1} - Y_t = (X_t)'G X_{t+1} .
```
Let $M_t = \exp(Y_t)$.
Define   a matrix ${\sf M}$ whose $(i,j)$th element is ${\sf M}_{ij} = \exp(G_{ij})`[^footnoteexp]. The stochastic process $\{M_t\}$ is governed by the *multiplicative functional*:
```{math}
:label: eqn:multfnMark
{\frac {M_{t+1}}{M_t}} = \exp\left[ (X_t)' G X_{t+1} \right] = (X_t)'{\sf M} X_{t+1}.
```
Associated with this multiplicative functional is the principal eigenvalue problem
```{math}
E\left[ {\frac {M_{t+1}}{M_t}}  \tilde e \cdot X_{t+1} | X_t = x \right] = \exp\left( \tilde \eta \right) \tilde e \cdot x.
```
Write the $j^{th}$ entry of ${\tilde e}$ as ${\tilde e}_j$.   Notice that because $X_t$ always assumes the value of one of the unit vectors $s_i, i =1, \ldots, n$, we can write
```{math}
(X_t)'{\sf M} X_{t+1} = {\sf M}_{ij}
```
when $X_t = s_i$ and $X_{t+1} = s_j$.
Therefore, we can write the principal eigenvalue problem as
```{math}
\sum_j P_{ij} {\sf M}_{ij} \tilde e_j = \exp(\tilde \eta) \tilde e_i
```
or
```{math}
:label: eqn:peigMarkov
\widetilde P \tilde e = \exp(\tilde \eta) \tilde e
```
where $\widetilde P_{ij} = P_{ij} {\sf M}_{ij}$.  We want the largest eigenvalue and  associated eigenvector of {eq}`eqn:peigMarkov`.
After solving the principal eigenvalue problem, compute
```{math}
:label: eqn:tildeMform
{\widetilde {\sf M}}_{ij} = \exp\left( - \tilde \eta \right) {\sf M}_{ij} \frac {e_j}{e_i}
```
and form ${\widetilde {\sf M}} = [{\widetilde {\sf M}}_{ij}]$.
By construction, ${\widetilde {\sf M}}$ is a stochastic  matrix that we can
use to form   increments $(X_t)'{\widetilde {\sf M}} X_{t+1}$ in a positive multiplicative martingale process $\{{\widetilde M_t}\}$:
```{math}
:label: eqn:likeMarkov
{\frac {{\widetilde M}_{t+1}}{{\widetilde M}_t}} = (X_t)'{\widetilde {\sf M}} X_{t+1}.
```
To achieve a {prf:ref}`prop:factor` representation of  the multiplicative functional $M_t$, use formula {eq}`eqn:tildeMform` for ${\widetilde{\sf M}}_{ij}$ to get
$
{\sf M}_{ij} = \exp\left( \tilde \eta \right) {\widetilde {\sf M}}_{ij}  \frac {e_i}{e_j} $
so that we can write {eq}`eqn:multfnMark` as
```{math}
:label: eqn:multplMarkfinal
{\frac {M_{t+1}}{M_t}} =  \exp\left( \tilde \eta \right) \left[(X_t)'{\widetilde {\sf M}} X_{t+1}\right]
\left(  \frac {e \cdot X_t }{e \cdot X_{t+1}} \right) .
```

[^footnoteexp]: This construction of
${\sf M}$ exploits the fact that $X_t$ is a coordinate vector.
````

````{prf:remark}
Let $i_t$ be the index of the Markov state at time $t$ and let $\{ i_0, i_1, \ldots, i_T\}$ be a simulation of the (asymptotically stationary) Markov process for $\{X_t\}$. We can use {eq}`eqn:multfnMark` to generate a simulation of the multiplicative functional $M_t$ and equation {eq}`eqn:likeMarkov` to generate a simulation of the positive multiplicative martingale $\widetilde M_t$.
````

*Construct some Python examples here and nearby*


````{prf:theorem}
:label: fact:MC_mult_forecasting
(Forecasting a Multiplicative Functional, II)

Let $\{M_t\}$ be the multiplicative functional described in {prf:ref}`ex:Markov_chain_mult` with transition matrix $P$ and matrix $\mathsf{M}$ defining multiplicative increments. Form a matrix $\widetilde P$ whose $i,j$ element is $P_{ij} \mathsf{M}_{ij}$. Then for integer $j \geq 1$
```{math}
:label: eqn:MC_forecast
E \left[ M_{t+j} | M_t, X_t = s_i \right] =    \left( \sum_k \widetilde P^{(j)}_{i,k} \right) M_t ,
```
where $\widetilde P^{(j)}$ is the $j$th power of $\widetilde P$.
````

````{prf:example}
:label: ex:Markchain_Lucas_1

(Stochastically growing payouts)

Assume the same Markov chain setting of example {prf:ref}`ex:Markov_chain_mult`.
Define a matrix $\sf{S}$ whose $(i,j)$th element is $\sf{S}_{ij} = \exp(G_{S,ij})$, where $G_{S,ij}$ is a stochastic discount rate for moving from state $s_i$ at time $t$ to state $s_j$ at time $t+1$.[^Sconstruction] A stochastic discount factor process $\{S_t\}_{t=0}^\infty$ is governed by the *multiplicative functional*:

```{math}
{\frac {S_{t+1}}{S_t}} = \exp\left[ (X_t)' G_S X_{t+1} \right] = (X_t)'{\sf S} X_{t+1}.
```

Define a matrix $\sf{D}$ whose $(i,j)$th element is $\sf{D}_{ij} = \exp(G_{D,ij})$.
A non-negative payout or dividend process $\{d_t\}_{t=0}^\infty$ is governed by the *multiplicative functional*:

```{math}
{\frac {d_{t+1}}{d_t}} = \exp\left[ (X_t)' G_D X_{t+1} \right] = (X_t)'{\sf D} X_{t+1}.
```

As in example {prf:ref}`ex:Markov_chain_mult`, because $X_t$ always assumes the value of one of the unit vectors $s_i, i =1, \ldots, n$, we can write

$(X_t)'{\sf S} X_{t+1} = {\sf S}_{ij}$   and $(X_t)'{\sf D} X_{t+1} = {\sf D}_{ij}$

when $X_t = s_i$ and $X_{t+1} = s_j$.

Let $p_t$ be the price at the beginning of period $t$ of a claim to the stochastically growing or shrinking stream of payouts $\{d_{t+j}\}_{j=0}^\infty$. It satisfies

```{math}
p_t = E\left[\frac{S_{t+1}}{S_t} (d_t + p_{t+1})\right] \Bigl| \mathfrak{F}_t ,
```

or

```{math}
\frac{p_t}{d_t} =
E\left[\frac{S_{t+1}}{S_t}
       \left(1 + \frac{d_{t+1}}{d_t} \frac{p_{t+1}}{d_{t+1}}\right)
 \right] \Bigl| \mathfrak{F}_t  .
```

where time $t$ information set $\mathfrak{F}_t$ includes $X_t, S_t, d_t$.
Guessing that the price-dividend ratio $\frac{p_t}{d_t}$ is a function of the Markov state $X_t$ only, and letting it equal $v_i$ when $X_t = s_i$, write the preceding equation as

```{math}
v_i =   \sum_{j=1}^n P_{ij} \left[  {\sf S}_{ij} \mathbf{1}  + \left( {\sf S}_{ij} \cdot {\sf D}_{ij} \right) v_j \right]
```

or

```{math}
v = c + \widetilde{P} v ,
```

where $c = \widehat{P} \mathbf{1}$ is by construction a nonnegative vector and we have defined nonnegative matrices $\widetilde{P} \in \mathbb{R}^{n \times n}$ and $\widehat{P} \in \mathbb{R}^{n \times n}$ by

```{math}
\begin{aligned}
\widetilde{P}_{ij} &= P_{ij} {\sf S}_{ij} {\sf D}_{ij}, \\
\widehat{P}_{ij} &= P_{ij} {\sf S}_{ij},
\end{aligned}
```

The equation $v = \widetilde{P} v + c$ has a nonnegative solution if and only if all eigenvalues of $\widetilde{P}$ are smaller than $1$ in modulus.
A sufficient condition for this is that all column sums, or that all the row sums, of $\widetilde{P}$ are less than one, which holds when $G_S + G_D \ll 0$. This condition asserts that discounting overtakes growth in dividends.
Given a solution $v$, the price-dividend ratio is a stationary process that is a fixed function of the Markov state:

```{math}
\frac{p_t}{d_t} = v_i \text{ when $X_t = s_i$}.
```

Meanwhile, the asset price process and the dividend process are both multiplicative functionals that exhibit either multiplicative growth or decay.

[^Sconstruction]: This construction of $\sf{S}$ again exploits the fact that $X_t$ is a coordinate vector.

````

````{prf:example}
:label: ex:Markchain_Lucas_2
(Lucas asset pricing model with growth)

Consider the setting of example {prf:ref}`ex:Markchain_Lucas_1`.
Let $d_t = C_t$, aggregate consumption, and let
```{math}
\frac{S_{t+1}}{S_t} = \exp(-\delta) \left(\frac{C_{t+1}}{C_t} \right)^{-\gamma} ,
```
where $\delta > 0$ is a rate of time preference and $\gamma \geq 1$ is a coefficient of relative risk aversion.
To obtain this special case of example {prf:ref}`ex:Markchain_Lucas_1`, we set
```{math}
\textsf{S}_{ij} = \exp(-\delta) \textsf{D}_{ij}^{-\gamma} ,
```
where we now interpret  $\textsf{D}_{ij}$ as the multiplicative  rate of growth of the level of aggregate consumption between $t$ and $t+1$
when $X_t = s_i$ and $X_{t+1} = s_j$.
````

````{prf:example}
:label: lab:term_struct_MC
(Term structure of interest rates, II)

Consider the framework of example {prf:ref}`ex:Markchain_Lucas_1` or {prf:ref}`ex:Markchain_Lucas_2`.
When the Markov state $X_t = s_i$ at time  $t$, the price of a risk-free zero-coupon bond paying one unit of consumption at time $t+j$
is
```{math}
p_{j,t} = E \left[ \frac{S_{t+j}}{S_t}  \right] \Big| X_t = s_i .
```
Let $\widehat{P}_{ij} = P_{ij} {\sf S}_{ij}$ and apply formula {eq}`eqn:MC_forecast` from fact {prf:ref}`fact:MC_mult_forecasting`
to deduce
```{math}
p_{j,t} =    \left( \sum_k \widehat P^{(j)}_{i,k} \right) ,
```
where $\widehat P^{(j)}$ is the $j$th power of the matrix $\widehat P$.
````

## Technical details involving multiplicity

When the [eigenfunction problem](prob:eigen) has multiple solutions, there is a unique solution
for which the process $\{ X_t \}$ is stationary and ergodic under the implied change of measure, namely,
 the solution associated with the minimum eigenvalue.

````{prf:theorem}
:label: prop:unique
Let $\{ M_t \}$ be a multiplicative functional. Suppose that $(\tilde \eta, \tilde e)$ solves [eigenfunction problem](prob:eigen) and that under the change of measure $\widetilde P$ implied by the associated martingale $\{ \widetilde M_t \}$ the stochastic process $\{ X_t \}$ is stationary and ergodic. Consider any other solution $(\eta^*, e^*)$ to [eigenfunction problem](prob:eigen) with implied martingale $\{ M_t^* \}$. Then

1. $\eta^* \ge \tilde \eta$.

2. If $\{ X_t \}$ is stationary and ergodic under the change of measure $Pr^*$ implied by the martingale $\{ M_t^* \}$, then $\eta^* = \tilde \eta$, $e^*$ is proportional to $\tilde e$, and $M_t^* = \widetilde M_t$ for all $t=0,1,... $.

````

First we show that $\eta^* \ge \tilde \eta$. Write:

```{math}
\mathbb{M}^t{e^*}(x) = \exp\left( \tilde \eta t \right) \widetilde{E} \left( \left[ \frac { \tilde e(X_0)}{\tilde e(X_t) } \right] {e^*}(X_t) \vert X_0=x \right) = \exp \left( \eta^* t \right) e^*(x) .
```

Thus,

```{math}
\widetilde{E} \left( \left[ \frac { e^*(X_t)}{\tilde e(X_t) } \right] \Biggl| X_0=x \right) = \exp \left( \eta^* t - \tilde \eta t \right) \left[ \frac {e^*(x)}{\tilde e(x)} \right] .
```

If $\tilde \eta > \eta^*$, then

```{math}
\lim_{t \rightarrow \infty} \widetilde{E} \left( \left[ \frac { e^*(X_t)}{\tilde e(X_t) } \right] \Biggl| X_0=x \right) = 0.
```

But this equality cannot be true because $\widetilde{Pr}$ implies that $\{ X_t \}$ is stationary and $\frac {e^*}{\tilde e}$ is strictly positive. Therefore, inequality must be satisfied.

Consider next the case in which $\eta^* > \tilde \eta$. Write

```{math}
\frac {M_t}{ M_0} = \exp\left( \eta^* t \right) \left( \frac { M^*_t}{ M^*_0 } \right) \left( \frac {e^*(X_0)}{e^*(X_t) } \right),
```

which implies that

```{math}
\mathbb{M}^t\tilde e(x) = \exp\left( \eta^* t \right) E^* \left[ \left( \frac {e^*(X_0)}{e^*(X_t) } \right) \tilde e(X_t) | X_0=x \right] = \exp \left( \tilde \eta t \right) \tilde e(x) .
```

Thus,

```{math}
E^* \left( \left[ \frac { \tilde e(X_t)}{e^*(X_t) } \right] \Biggl| X_0=x \right) = \exp\left( \tilde \eta t - \eta^* t \right)\left[ \frac {\tilde e(x)}{e^*(x)}\right].
```

Suppose that $\tilde \eta < \eta^*$, then

```{math}
\lim_{t \rightarrow \infty} E^* \left( \left[ \frac { \tilde e(X_t)}{e^*(X_t) } \right] \Biggl| X_0=x \right) = 0 ,
```

so that $\{ X_t : t=0,1,... \}$ cannot be stationary under the $Pr^*$ measure. Next suppose that $\tilde \eta = \eta^*$ and that $\frac {\tilde e(x)}{e^*(x)}$ is not constant. Then

```{math}
E^* \left( \left[ \frac { \tilde e(X_t)}{e^*(X_t) } \right] \Biggl| X_0=x \right) = \frac {\tilde e(x)}{e^*(x)}
```

and $\{ X_t : t=0,1,... \}$ cannot be ergodic under the $Pr^*$ measure.
## Long-run Growth Rates

It is instructive to use the Proposition {prf:ref}`prop:factor` factorization to study the behavior of ${\mathbb M}^jf$ for large $j$. Suppose that

```{math}
:label: eqn:fcondition
f > 0 \textrm{ and } 0 < {\widetilde E}\left[{\frac {f(X_t)} {{\tilde e}(X_t)}}\right] < \infty.
```

Then if $\{ X_t\}$ is stationary and ergodic under ${\widetilde P}r$

```{math}
{\frac 1 j} \log {\mathbb M}^jf(x)  =
{\tilde \eta} + {\frac 1 j} \log {\widetilde E} \left[  {\frac {f(X_j)} {{\tilde e}(X_j)}} \Big| X_0 = x \right]  - {\frac 1 j} \log {\tilde e} (x)
 \le  {\tilde \eta}  + {\frac 1 j} \log \sum_{t=1}^j {\widetilde E} \left[  {\frac {f(X_t)} {{\tilde e}(X_t)}} \Big| X_0 = x \right]
 - {\frac 1 j} \log {\tilde e}(x) \le {\tilde \eta}  +  {\frac 1 j} \log j + {\frac 1 j} \log {\frac 1 j} \sum_{t=1}^j
{\widetilde E} \left[  {\frac {f(X_t)} {{\tilde e}(X_t)}} \Big| X_0 = x \right] - {\frac 1 j} \log {\tilde e}(x) .
```

Consequently,

```{math}
\limsup_{j \rightarrow \infty} {\frac 1 j} \log {\mathbb M}^jf(x) \le {\tilde \eta} ,
```

so ${\tilde \eta}$ is an upper bound on the growth rate for a class of functions $f$ satisfying condition {eq}`eqn:fcondition`. The bound is attained when $f = {\tilde e}$. If in addition,

```{math}
:label: liminf
\liminf_{j \rightarrow \infty}  {\widetilde E}\left[{\frac {f(X_j)} {{\tilde e}(X_j)}} \Big| X_0 = x \right] > 0,
```

then

```{math}
\liminf_{j \rightarrow \infty}   {\frac 1 j} \log {\mathbb M}^jf(x) = {\tilde \eta}
```

and therefore

```{math}
\lim_{j \rightarrow \infty} {\frac 1 j} \log {\mathbb M}^jf(x) = {\tilde \eta}.
```

We now consider the consequences of strengthening the ergodicity restriction with a stochastic stability condition.


````{prf:condition}
:label: stochasticstability

A process $\{ X_t \}$ is **stochastically stable** under ${\widetilde {Pr}}$ if it is stationary for almost all $x$ under the distribution for $X_0$ induced by ${\widetilde {Pr}}$ 
```{math}
\lim_{j \rightarrow  \infty} {\widetilde E} \left[h(X_j) \mid X_0 = x \right] = {\widetilde E} \left[ h(X_t) \right]
```
for any Borel measurable $h$ satisfying ${\widetilde E} |h(X_t)| < \infty$.
````

Under condition {prf:condition}`stochasticstability`, we obtain a more refined approximation:

```{math}
\lim_{j \rightarrow \infty} \exp (- {\tilde \eta} j ) {\mathbb M}^j f(x)  = \lim_{j \rightarrow \infty}  {\widetilde E}
\left[  {\frac {f(X_j)} {{\tilde e}(X_j)}}  \Big| X_0 = x \right] {\tilde e}(x)
 =  {\widetilde E}\left[{\frac {f(X_t)} {{\tilde e}(X_t)}}\right] \tilde e(x).
```

Once we adjust for the impact of ${\tilde \eta}$, the limiting function is proportional to ${\tilde e}$. The function $f$ determines only a scale factor ${\widetilde E}\left[{\frac {f(X_t)} {{\tilde e}(X_t)}}\right] \tilde e(x)$.
## Applications and Extensions

We will apply and extend results in the chapter in a variety of ways in chapters [Likelihood-Based Inference](like) and [Elasticity and its Applications](elastic).

[^TSNote]: Lay out road map here.

