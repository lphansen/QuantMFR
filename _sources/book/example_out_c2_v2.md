(chap:markov)=
# Markov Processes

We call a random vector $X_t$ the *state* because it  describes probabilistically the position of a dynamic system at time $t$ from the perspective of a model builder or an econometrician. We construct a consistent sequence of probability distributions $Pr_\ell$ for a sequence of random vectors

```{math}
X^{[\ell]} \doteq \begin{bmatrix} X_0 \\ X_1 \\ \vdots \\ X_\ell \end{bmatrix}
```

for all nonnegative integers $\ell$ by specifying the following two elementary components of a *Markov process*: (i) a probability distribution for $X_0$, and (ii) a time-invariant distribution for $X_{t+1}$ conditional on $X_t$ for $t \geq 0$.   The vector $X_t$ suffices for conditioning on the past history of the process.  All other probabilities are functions of these two distributions. By creatively defining the state vector $X_t$, a Markov specification includes many models used in applied research.
## Constituents

Assume a state space $\mathcal{X}$ and a transition distribution $P(dx^*|x)$. For example, $\mathcal{X}$ could be $\mathbb{R}^n$ or a subset of $\mathbb{R}^n$.
The transition distribution $P$ is a conditional probability measure for each $X_t = x$ in the state space,
so it satisfies

```{math}
\int_{\{x^* \in \mathcal{X} \}}P(dx^* | x) = 1
```


for every $x$ in the state space.
If in addition we specify a marginal distribution $Q_0$ for the initial state $x_0$ over $\mathcal{X}$, then we have completely specified all joint distributions for the
stochastic process $\{X_t, t = 0, 1, \ldots\}$.  
The notation $P(dx^*|x)$ denotes a conditional probability measure; integration is over $x^*$ and conditioning is captured by $x$.
Thus, $x^*$ is a possible realization of next period's state and $x$ is a realization of this period's state.
The conditional probability measure $P(dx^* |x)$ assigns conditional probabilities to next period's state given that this period's state is $x$.
Often, but not always, the conditional distributions have densities against a common distribution $\lambda(dx^*)$ to be used to integrate over states.
That lets us use a *transition density* to represent the conditional probability measure.
````{prf:example}
:label: ex:ex01
A first-order vector autoregression is a Markov process. In this example, we consider $n$ such processes, indexed by $i$. 
The index $i$ represents a discrete form of parameter uncertainty, and we may include $i$ as an additional  time-invariant component of the state.  
Here  $P(dx^*|x, i )$ is a normal distribution with mean ${\mathbb A}_ix$ and covariance matrix ${\mathbb B}_i{{\mathbb B}_i}'$ for square matrices ${\mathbb A}_i$ and  matrices ${\mathbb B}_i$ with full column rank.[^singularity] These assumptions imply the vector autoregressive (VAR) representation for $i=1,2,..., n.$
```{math}
X_{t+1} = {\mathbb A}_i X_t + {\mathbb B}_i W_{t+1} ,
```
for $t \geq 0$, where $W_{t+1}$ is a multivariate standard normally distributed random vector that is independent of $X_t, i$.

[^singularity]: When ${\mathbb B}_i{\mathbb B}_i'$ is singular, a density may not exist with respect to Lebesgue measure. The covariance matrix ${\mathbb B}_i{\mathbb B}_i'$ is typically singular for a first-order vector autoregression constructed by rewriting a higher-order vector autoregression.
````
````{prf:example} 
:label: ex:ex02 
A discrete-state Markov chain consists of a 
$Q_0$ represented as a row vector and a transition probability $P(dx^*|x)$  represented as a matrix  with one row and one column for each possible value of the state $x$. 
Rows contain vectors of probabilities of next period's state conditioned on a realized value of 
this period's state.  As with the VAR example, we include $n$ such process and augment the state with a time-invariant component that captures which transition matrix captures the actual transitional dynamics.   
````

It is useful to construct an operator by applying a one-step conditional expectation operator to functions of a Markov state. Let $f:{\mathcal X} \rightarrow {\mathbb R}$.
For bounded $f$, define:

```{math}
:label: eqn:Toperatordef
{\mathbb T} f (x) = E \left[ f(X_{t+1}) | X_t = x \right] = \int_{\{x^* \in {\mathcal X}\}}  f(x^*) P(d x^*|x).
```

The Law of Iterated Expectations justifies iterating on ${\mathbb T}$  to form conditional expectations of the function $f$ of the Markov state over longer horizons:

```{math}
{\mathbb T}^j f(x) = E  \left[ f(X_{t+j}) | X_t = x \right].
```

The operator ${\mathbb T}$  gives an alternative way to represent the transitional dynamics for a Markov process. Indeed, by applying ${\mathbb T}$ to a suitable range of test functions $f$, we can construct a conditional probability measure.  

````{prf:theorem}
Start with a conditional expectation operator ${\mathbb T}$ that maps a space of bounded functions into itself. We can use ${\mathbb T}$ to construct a conditional probability measure $P(dx^*|x)$ provided that ${\mathbb T}$ is (a) well-defined on the space of bounded functions, (b) preserves the bound, (c) maps nonnegative functions into nonnegative functions, and (d) maps the unit function into the unit function.
````
## Stationarity

We can construct a stationary Markov process by carefully choosing the distribution of the initial state $X_0$.
````{prf:definition}
:label: def:stationdist

A probability measure $Q$ over a state space $\mathcal{X}$ for a Markov process with transition probability $P$ is a **stationary distribution** if it satisfies

```{math}
\int_{ \{ x \in \mathcal{X} \}} P(d x^*|x) Q(dx) = Q(d x^*).
```
````

We will sometimes refer to a stationary density $q$. A density is always relative to a measure. With this in mind, let $\lambda$ be a measure used to integrate over possible Markov states on the state space $\mathcal{X}$. Then a density $q$ is a nonnegative (Borel measurable) function of the state for which $\int q(x) \lambda(dx) = 1$.
````{prf:definition}
A **stationary density** over a state space $\mathcal{X}$ for a Markov process 
with transition probability $P$ is a probability density $q$ with respect to a measure $\lambda$ over the state space $\mathcal{X}$ that satisfies
```{math}
\int P(d x^*|x) q(x) \lambda(dx) = q(x^*) \lambda(dx^*).
```
````
````{prf:definition}
:label: def:reversible

A Markov process with stationary density $q$ and transition density $P(dx^*|x)$ is said to be **reversible** if

```{math}
:label: eqn:reversible
P(dx^*|x)q(x) \lambda(d x) = P(dx|x^*)q(x^*)  \lambda(d x^*).
```
````
````{prf:example}
Various sufficient conditions imply the existence of a stationary distribution. Given a transition distribution $P$, one such condition that is widely used to justify some calculations from numerical simulations is that the Markov process be *time reversible*, which means that

```{math}
:label: eqn:reversiblenew
P(dx^*|x) Q(dx) = P(dx|x^*) Q(dx^*)
```

for some probability distribution $Q$ on $\mathcal{X}$. Because a transition distribution satisfies $\int_{\{ x \in \mathcal{X}\}}  P(dx|x^*) =1 $,

```{math}
\int_{\{ x \in \mathcal{X}\}} P(dx^*|x) Q(dx)  = \int_{\{ x \in \mathcal{X}\}} P(dx|x^*) Q(dx^*)  = Q(dx^*) ,
```
so $Q$ is a stationary distribution by {prf:ref}`def:stationdist`. Restriction {eq}`eqn:reversiblenew` implies that the process is time reversible in the sense that forward and backward transition distributions coincide. Time reversibility is special, so later we will explore other sufficient conditions for the existence of stationary distributions.[^bayesianMarkov]

````

[^bayesianMarkov]: Numerical Bayesian statistical analysis often computes a posterior probability distribution by iterating to convergence a reversible Markov process whose stationary distribution is that posterior distribution.
````{prf:remark}
When a Markov process starts at a stationary distribution, we can construct the process $\{ X_t : t=1,2,...\}$ with a measure-preserving transformation ${\mathbb S}$ of the type featured in chapter [](chap:process), section [](sec:stochprocessconstructionI).
````

%\textcolor{red}{Tom XXXXX: add an example here.}
%
Given a stationary distribution $Q$, form the space of functions ${\mathcal L}^2$

```{math}
{\mathcal L}^2 = \{ f:{\mathcal X} \rightarrow {\mathbb R} : \int f(x)^2 Q(dx)  < \infty \} .
```

It can be shown that ${\mathbb T} : {\mathcal L}^2 \rightarrow {\mathcal L}^2$. On this space, a well-defined norm is

```{math}
\| f \| = \left[\int f(x)^2 Q(dx)\right]^{1/2} .
```
(sec:eigfns)=
## ${\mathcal L}^2$ and Eigenfunctions

We connected ergodicity to a statistical notion of invariance in chapter [](chap:process).  The word invariance brings to mind a generalization of eigenvectors called eigenfunctions. Eigenfunctions of a linear mapping characterize an invariant subspace of functions such that the application of a linear mapping to any element of that space remains in the same subspace.  Eigenfunctions associated with a unit eigenvalue are themselves invariant under the mapping.  So perhaps it is not surprising that such eigenfunctions of ${\mathbb T}$ come in handy for studying ergodicity of Markov processes.  

Given a stationary distribution $Q$, form the space of functions 
```{math}
{\mathcal L}^2 = \{ f:{\mathcal X} \rightarrow {\mathbb R} : \int f(x)^2 Q(dx)  < \infty \} .
```
It can be verified that ${\mathbb T} : {\mathcal L}^2 \rightarrow {\mathcal L}^2$ and that 
```{math}
\| f \| = \left[\int f(x)^2 Q(dx)\right]^{1/2} 
```
is a well-defined norm on ${\mathcal L}^2$.

We now study eigenfunctions of the conditional expectation operator ${\mathbb T}$.
````{prf:definition}
A function $f \in \mathcal{L}^2$ that solves  $\mathbb{T}  f =  f$
is  an eigenfunction of $\mathbb{T}$  associated with a unit eigenvalue.[^eigenfunctionNote]

[^eigenfunctionNote]: For Markov processes, all invariant events depend only on the initial $X$.  A reference  is {cite:t}`doob`, Theorem 1.1, page 460.  Indicator functions of these are events are thus representable as eigenfunctions associated with unit eigenvalues.
````

The following proposition asserts that an eigenfunction  $\tilde{f}(X_t)$ associated with a unit eigenvalue is constant as $X_t$ moves through time.
````{prf:proposition}
:label: lem:uniteigen
Suppose that $\tilde{f}$ is an eigenfunction of $\mathbb{T}$ associated with a unit eigenvalue. Then $\{\tilde{f}(X_t) : t=0,1,...\}$ is constant over time with probability one.
````


````{prf:proof}


```{math}
E \left[\tilde{f}(X_{t+1}) \tilde{f}(X_t)\right] = \int (\mathbb{T}\tilde{f})(x) \tilde{f}(x) Q(dx) = \int \tilde{f}(x)^2 Q(dx) =
E \left[\tilde{f}(X_t)^2\right]
```
where the first equality follows from the Law of Iterated Expectations. Then because $Q$ is a stationary distribution,
```{math}
\begin{eqnarray*}
E\left([\tilde{f}(X_{t+1}) - \tilde{f}(X_t)]^2\right)  & = & E\left[\tilde{f}(X_{t+1})^2\right] + E \left[\tilde{f}(X_t)^2\right] \cr &&- 2 E\left[ \tilde{f}(X_{t+1})\tilde{f}(X_t) \right] \cr
  & = & 0. \end{eqnarray*}
```
````

(sec:MarkErgodic)=
## Ergodic Markov Processes

Chapter [](chap:process) studied special statistical models that, because they are ergodic, are affiliated with a Law of Large Numbers in which limit points are constant across sample points $\omega \in \Omega$. Section [](sec:ergodic_decomp) described other statistical models that are not ergodic and that are components of more general probability specifications that we used to express the idea that a statistical model is unknown.[^unknownParameters] As we described, even when the statistical model is unknown, ergodic processes remain of interest as they are building blocks (specific statistical models) that are revealed over time.  We now explore ergodicity in the context of Markov processes.

From {prf:ref}`lem:uniteigen` we know that time-series averages of an eigenfunction ${\mathbb T} \tilde f = \tilde f$ are invariant over time, so

```{math}
{\frac 1 N} \sum_{t=1}^N \tilde f(X_t) = \tilde f(X).
```

However, when ${\tilde f}(x)$ varies across sets of states $x$ that occur with positive probability under $Q$, a time series average ${\frac 1 N} \sum_{t=1}^N \tilde f(X_t)$ can differ from $\int \tilde f(x) Q(dx)$. This happens when observations of $\tilde f(X_t)$ along a sample path for $\{X_t\}$ convey an inaccurate impression of how $f(X)$ varies across the stationary distribution $Q(dx)$. See {prf:ref}`ex:MC2` below. We can exclude the possibility of such inaccurate impressions by imposing a restriction on the eigenfunction equation ${\mathbb T}f = f$. % to state a sufficient condition for ergodicity.


[^unknownParameters]: Unknown parameters manifest themselves as unknown statistical models.


````{prf:proposition}
:label: prop:ergo
When a unique solution to the equation
```{math}
{\mathbb T}f = f
```
is a constant function (with $Q$ measure one), then it is possible to construct $\{ X_t : t=0,1,2,...\}$ 
as a stationary and ergodic Markov process with ${\mathbb T}$ as 
the one-period conditional expectation operator and $Q$ as the initial distribution for
$X_0$.[^ergo_footnote]  
````

[^ergo_footnote]: In particular, the process can be represented using a probability measure $Pr$ defined over events in ${\mathfrak F}$, a transformation ${\mathbb S}$ for which$({\mathbb S}, Pr)$ is measure preserving, and ergodic and a measurement function ${\widetilde X}$ such that $\left\{ {\widetilde X}\circ {\mathbb S}^t : t=0,1, \ldots \right\}$ has the same induced distribution as the process $\{X_t : t=0,1,2, \ldots \}$.

Evidently, ergodicity is a property that obtains relative to a stationary distribution $Q$ of the Markov process. If there are multiple stationary distributions, it is possible that there is a unique constant function $f$ that solves ${\mathbb T}f = f$ problem for one stationary distribution and that non-constant solutions exist for other stationary distributions. 

### Invariant events for a Markov process

Consider an eigenfunction ${\tilde f}$ of ${\mathbb T}$ associated with a unit eigenvalue. Let $\phi : {\mathbb R} \rightarrow {\mathbb R}$ be a bounded Borel measurable function. Since $\{ {\tilde f}(X_t) : t=0,1,2,... \}$ is invariant over time, so is $\left\{ \phi\left[{\tilde f}(X_t)\right] : t=0,1,2, \ldots \right\}$ and it is necessarily true that

```{math}
{\mathbb T} (\phi \circ {\tilde f}) = \phi \circ {\tilde f}.
```
Therefore, from an eigenfunction ${\tilde f}$ associated with a unit eigenvalue, we can construct other eigenfunctions,[^eigen_footnote] for example
```{math}
:label: newjunk1
\phi[{\tilde f}(x)] = \begin{cases} 
1 & \text{if} \ {\tilde f}(x) \in {\tilde {\mathfrak b}} \\
0 & \text{if} \ {\tilde f}(x) \notin {\tilde {\mathfrak b}} \end{cases}
```
for some Borel set ${\tilde {\mathfrak b}}$ in ${\mathbb R}$.

[^eigen_footnote]: This construction also works for unbounded functions $\phi$ provided that $\phi \circ \tilde f$ is square integrable under the $Q$ measure.

It follows that
```{math}
\Lambda = \{ \omega \in \Omega : {\tilde f}[X_0(\omega)] \in {\tilde {\mathfrak b}} \}
```
is an invariant event in $\Omega$. Note that by constructing the Borel set, ${\mathfrak b}$ in $\mathcal X$
```{math}
{\mathfrak b} = \{ x : {\tilde f}(x) \in {\tilde {\mathfrak b}}  \}
```
we can represent $\Lambda$ as
```{math}
:label: invariantrep
\Lambda = \{ \omega \in \Omega : X_0(\omega) \in {\mathfrak b} \}.
```
Thus we have shown how to construct many non-degenerate eigenfunctions, starting from an initial such function.

For Markov processes, all invariant events can be represented as in {eq}`invariantrep`, which is expressed in terms of the initial state $X_0$. See {cite:t}`doob`. Thus, associated with an invariant event is a Borel set in ${\mathcal X}$. Let ${\mathfrak J}$ denote the collection of Borel subsets of ${\mathcal X}$ for which $\Lambda$ constructed as in {eq}`invariantrep` is an invariant event. From these invariant events, we can also construct many non-degenerate eigenfunctions as indicator functions of sets in ${\mathfrak J}$. Formally, if ${\mathfrak b} \in {\mathfrak J}$, then the indicator function
```{math}
:label: neweigen
f(x) = \begin{cases} 
1 & \text{if} \ x \in {\mathfrak b} \\
0 & \text{if} \ x \notin {\mathfrak b} \end{cases}
```
satisfies
```{math}
{\mathbb T} f = f
```
with $Q$ probability one. Provided that the probability of $\Lambda$ is neither zero nor one, then we have constructed a nonnegative function $f$ that is strictly positive on a set of positive $Q$ measure and zero on a set with strictly positive $Q$ measure.


More generally, when a Markov process $\left\{X_t: t \geq 0\right\}$ is not ergodic, there exist bounded eigenfunctions with unit eigenvalues that are not constant with $Q$ measure one. For a non-degenerate eigenfunction $\tilde{f}$ with unit eigenvalue to be constant with $Q$ measure one, it shouldn't be possible for the Markov process permanently to get stuck in a subset of the state space which has probability different from one or zero. 

Suppose now we consider any Borel set $\mathfrak{b}$ of $\mathcal{X}$ that has $Q$ measure that is neither zero nor one. Let $f$ be constructed as in {eq}`neweigen` without restricting $\mathfrak{b}$ to be in $\mathfrak{J}$. Then $\mathbb{T}^j$ applied to $f$ is the conditional probability of $\left\{X_j \in \mathfrak{b}\right\}$ as of date zero. If we want time series averages to converge to unconditional expectations, we must require that the set $\mathfrak{b}$ be visited eventually with positive probability. To account properly for all possible future dates we use a mathematically convenient resolvent operator defined by

```{math}
\mathbb{M} f(x)=(1-\lambda) \sum_{j=0}^{\infty} \lambda^j \mathbb{T}^j f.
```

for some constant discount factor $0<\lambda<1$. Notice that If $\tilde{f}$ is an eigenfunction of $\mathbb{T}$ associated with a unit eigenvalue, then the same is true for $\mathbb{T}^j$ and hence for $\mathbb{M}$. We translate the requirement that $\mathfrak{b}$ be eventually visited to a restriction that applying $\mathbb{M}$ the indicator function $f$ yields a strictly positive function. The following statement extends this restriction to all nonnegative functions that are distinct from zero.

````{prf:proposition}
:label: prop:ergod100
Suppose that for any $f \ge 0$
such that $\int f(x) Q(dx) > 0$, ${\mathbb M} f(x) > 0$ for all $x \in {\mathcal X}$ with $Q$ measure one.  Then
any solution ${\tilde f}$ to
${\mathbb T}  f =  f$ is necessarily constant with  $Q$ measure one. 
````

````{prf:proof}
Consider an eigenfunction ${\tilde f}$ associated with a unit eigenvalue.  The function $f = \phi \circ {\tilde f}$ necessarily satisfies:
```{math}
{\mathbb M}f = f
```
for any $\phi$ of the form {eq}`newjunk1`.  If such an $f$ also satisfies $\int f(x) Q(dx) > 0$, then $f(x)=1$
  with $Q$ probability one.  Since this holds for any Borel set ${\mathfrak b}$ in ${\mathbb R}$, ${\tilde f}$ must be constant with $Q$ probability one. 
````

{prf:ref}`prop:ergod100` supplies a sufficient condition for ergodicity. A more restrictive  sufficient condition is that there exists an integer $m \geq 1$ such that
 ```{math}
{\mathbb T}^{m} f(x) > 0
```
 for any $f \ge 0$ such that $\int f(x) Q(dx) > 0$
on a set with $Q$ measure one.

````{prf:remark}
The sufficient conditions imposed in {prf:ref}`prop:ergod100` imply a property called *irreducibility* relative to the probability measure $Q$. While this proposition presumes that $Q$ is a stationary distribution, *irreducibility* allows for a more general specification of $Q$.
````

{prf:ref}`prop:ergod100` provides a way to verify ergodicity. As discussed in chapter [](chap:process), ergodicity is a property of a statistical model. As statisticians or econometricians we often entertain a set of Markov models, each of which is ergodic. For each model, we can build a probability $Pr$ using the canonical construction given at the outset of chapter [](chap:process). Convex combinations of these probabilities are measure-preserving but not necessarily ergodic when used in conjunction with the shift transformation ${\mathbb S}$. We can take the ergodic Markov models to be the building blocks for a specification to be used in a statistical investigation. There can be a finite number of these building blocks or even a continuum of them represented in terms of an unknown parameter vector.
````{prf:definition}
The process $\{ X_t \}$ is said to be irreducible with respect to $\widetilde{Q}$ if for any $f \ge 0$ such that $\int f(x) \widetilde{Q}(dx) > 0$,
${\mathbb{M}} f(x) > 0$ for all $x \in {\mathcal{X}}$ with ${\widetilde{Q}}$ measure one.
````
````{prf:proposition}
When ${\widetilde Q}$ is a stationary distribution and $\left\{ X_t \right\}$ is irreducible with respect to $\widetilde Q$, the process is necessarily ergodic.
````

````{prf:proof}
By imitating the proof of {prf:ref}`prop:ergod100`, we can establish that irreducibility rules out bounded eigenfunctions that are not constant with ${\widetilde Q}$ measure one.
````
## Periodicity

Next, we study a notion of periodicity of a stationary and ergodic Markov process.[^periodicity_footnote] To define periodicity of a Markov process, for a given positive integer $p$ we construct a new Markov process by sampling an original process every $p$ time periods. This is sometimes called 'skip-sampling' at sampling interval $p$.[^skip_sampling_footnote] With a view toward applying {prf:ref}`lem:uniteigen` to ${\mathbb T}^p$, solve

```{math}
:label: perioddef
{\mathbb T}^p f = f
```

for a function ${\tilde f}$. We know from {prf:ref}`lem:uniteigen` that for an $\tilde f$ that solves {eq}`perioddef`, $\{ {\tilde f}(X_t) : t=0, p, 2p, \ldots \}$
is invariant and so is $\{ {\tilde f}(X_t) : t=1,p+1,2p+1,...\}$. The process ${\tilde f}(X_t)$ is periodic with period $p$ or $np$ for any positive integer $n$.

[^periodicity_footnote]: Our definition of periodicity is confined to a stationary distribution. Actually, periodicity can be defined more generally. We limit our treatment of periodicity to specifications of transition probabilities for which there exist stationary distributions for convenience here.

[^skip_sampling_footnote]: See {cite:t}`Hansen_Sargent_1993` and {cite:t}`HansenSargent_Recursive_Models`.
````{prf:definition}
The *periodicity* of an irreducible Markov process $\left\{ X_t \right\}$ with respect to ${\widetilde Q}$ is the smallest positive integer $p$ such that there is a solution to equation {eq}`perioddef` that is not constant with ${\widetilde Q}$ measure one. When there is no such integer $p$, we say that the process is *aperiodic*.
````

````{prf:result}
:label: result:pinterval
Consider a counterpart of the resolvent operator ${\mathbb M}$ constructed by sampling at interval given by positive integer $p$:

```{math}
:label: eqn:Mpdef
{\mathbb M}_p f(x) = (1 - \lambda) \sum_{j=0}^\infty \lambda^{j} {\mathbb T}^{pj} f.
```

Provided that ${\mathbb M}_p f(x) > 0$ with ${Q}$ measure one and all $p \ge 0$ for any $f \ge 0$ such that $\int f(x) Q(dx) > 0$, the Markov process is aperiodic.
````



[^SufficientConditionNote]: A more restrictive sufficient condition is again that there exists an $m$ such that ${\mathbb T}^{m} f(x) > 0$ on a set of $Q$ measure one for any $f \ge 0$ such that $\int f(x) Q(dx) > 0$. Given that this property holds for ${\mathbb T}^{m}$, it must also hold true for $p m$ for any nonnegative integer $p$.
(subsec:chain)=
## Finite-State Markov Chains

Suppose that $\mathcal{X}$ consists of $n$ possible states. We can label these states in a variety of ways, but for now we suppose that state $x_j$ is the coordinate vector consisting entirely of zeros except in position $j$, where there is a one. Let $\mathbb{P}$ be an $n$ by $n$ transition matrix, where entry ${i,j}$ is the probability of moving from state $i$ to state $j$ in a single period. Thus, the entries of $\mathbb{P}$ are all nonnegative and 
```{math}
\mathbb{P} \textbf{1}_n = \textbf{1}_n,
```
where $\textbf{1}_n$ is an $n$-dimensional vector of ones.

Let $\textbf{q}$ be an $n$-dimensional vector of probabilities. Stationarity requires that
```{math}
:label: qstab
\textbf{q}'\mathbb{P} = \textbf{q}',
```
where $\textbf{q}$ is a row eigenvector (also called a left eigenvector) of $\mathbb{P}$ associated with a unit eigenvalue.

We use a vector $\textbf{f}$ to represent a function from the state space to the real line. Each coordinate of $\textbf{f}$ gives the value of the function at the corresponding coordinate vector. Then the conditional expectation operator $\mathbb{T}$ can be represented in terms of the transition matrix $\mathbb{P}$:
```{math}
E(\textbf{f} \cdot X_{t+1} | X_t = x) = (\mathbb{T}\textbf{f}) \cdot x = x'\mathbb{P} \textbf{f}.
```
Now consider column eigenvectors called "right eigenvectors" of $\mathbb{P}$ that are associated with a unit eigenvalue.
````{prf:theorem}
:label: prop:finiteP1

Suppose that the only solutions to

```{math}
{\mathbb T} {\bf f} = {\bf f}
```
are of the form ${\bf f} \propto \textbf{1}_n$, where $\propto$ means 'proportional to'.
Then we can construct a process that is stationary and ergodic by initializing the process with density ${\bf q}$ determined by equation {eq}`qstab`. 
````

We can weaken the {prf:ref}`prop:finiteP1` sufficient condition for stationarity and ergodicity to allow nonconstant right eigenvectors. This weakening is of interest when there are multiple stationary distributions.
````{prf:theorem}
:label: prop:finiteP2
  Assume that there exists a real number $\mathbf{r}$ such that the right eigenvector $\mathbf{f}$ and a stationary distribution $\mathbf{q}$ satisfy
```{math}
\min_{\sf r} \sum_{i=1}^n ({\sf f}_i - {\sf r})^2 {\sf q}_i = 0.
```
Then the process is stationary and ergodic.
````

Notice that if ${\sf q}_i$ is zero, the contribution of ${\sf f}_i$ to the least squares objective can be neglected. This allows for non-constant $\mathbf{f}$'s, albeit in a limited way.

Three examples illustrate ideas in these propositions.
````{prf:example}
:label: ex:MC1

Recast {prf:ref}`ex:period` as a Markov chain with transition matrix
${\mathbb P}=\begin{bmatrix}0 & 1 \cr 1 & 0\end{bmatrix}$. This chain has a unique stationary distribution $ q=\begin{bmatrix}.5 & .5 \end{bmatrix}'$ and the invariant functions are $\begin{bmatrix} {\sf r} & {\sf r} \end{bmatrix}'$ for any scalar ${\sf r}$. Therefore, the process initiated from the stationary distribution is ergodic. The process is periodic with period two since the matrix ${\mathbb P}^2$ is an identity matrix and all two dimensional vectors are eigenvectors associated with a unit eigenvalue.
````

%{\bf {XXXXX Lars added the last two sentences. Perhaps we should switch orders of the examples.}}
````{prf:example}
:label: ex:MC2

Recast {prf:ref}`ex:invariant` as a Markov chain with transition matrix
${\mathbb P}=\begin{pmatrix}1 & 0 \\ 0 & 1\end{pmatrix}$. This chain has a continuum of stationary distributions $\pi \begin{pmatrix}1 \\ 0 \end{pmatrix}+ (1- \pi )\begin{pmatrix}0 \\ 1 \end{pmatrix}$ for any $\pi \in [0,1]$ and invariant functions $\begin{pmatrix} {\sf r}_1 \\ {\sf r}_2 \end{pmatrix}$ for any scalars ${\sf r}_1, {\sf r}_2$. Therefore, when $\pi \in (0,1)$ the process is not ergodic because if ${\sf r}_1 \ne {\sf r}_2$ the resulting invariant function fails to be constant across states that have positive probability under the stationary distribution associated with $\pi \in (0,1)$. When $\pi \in (0,1)$, nature chooses state $i=1$ or $i=2$ with probabilities $\pi, 1-\pi$, respectively, at time $0$. Thereafter, the chain remains stuck in the realized time $0$ state. Its failure ever to visit the unrealized state prevents the sample average from converging to the population mean of an arbitrary function of the state.
````
````{prf:example}
:label: ex:MC3

A Markov chain with transition matrix
```{math}
{\mathbb P}=\begin{bmatrix}.8 & .2 & 0  \cr .1  & .9 & 0 \cr
               0 & 0 & 1\end{bmatrix}
```
has a continuum of stationary distributions
```{math}
\pi \begin{bmatrix} {1\over 3} & {2 \over 3} & 0 \end{bmatrix}'
+(1- \pi) \begin{bmatrix} 0 & 0 & 1 \end{bmatrix}'
```
for $\pi \in [0,1]$ and invariant functions
```{math}
\begin{bmatrix} {\sf r}_1  &  {\sf r}_1 & {\sf r}_2 \end{bmatrix}'
```
for any scalars ${\sf r}_1, {\sf r}_2$. Under any stationary distribution associated with $\pi \in (0,1)$,
the chain is not ergodic because some invariant functions are not constant with probability one. But under stationary distributions associated with $\pi =1$ or $\pi=0$, the
chain is ergodic.
````
## Limited Dependence

Recall the conditional expectations operator ${\mathbb T}$ defined in equation {eq}`eqn:Toperatordef` for a space ${\mathcal L}^2$ of functions $f$ of a Markov process with transition probability $P$ and stationary distribution $Q$ and for which $f(X_t)$ has a finite second moment under $Q$:
```{math}
{\mathbb T} f (x) = E \left[ f(X_{t+1}) \mid  X_t = x \right] = \int_{\{x^* \in {\mathcal X}\}}  f(x^*) P(d x^*|x) .
```
We suppose that under the stationary distribution $Q$, the process is ergodic.  

Because it is often useful to work with random variables that have been 'centered' by subtracting out their means, we define the following subspace of ${\mathcal L}^2$:
```{math}
:label: def:N
{\mathcal N} = \left\{ f \in {\mathcal L}^2 :  \int f(x) Q(dx)   = 0 \right\}.
```
We use the same norm $\| f \| = \left[ \int f(x)^2 Q(dx)\right]^{1/2}$ on both ${\mathcal L}^2$ and ${\mathcal N}$ too.
````{prf:theorem} 
The conditional expectation operator $\mathbb{T}$ is said to be a *strong contraction* on $\mathcal{N}$ if there exists $0 < \rho < 1$ such that
```{math}
\| \mathbb{T} f \| \le \rho \| f \|
```
for all $f \in \mathcal{N}$.
````

When $\mathbb{T}^m$ is a strong contraction for some positive integer $m$ and some $\rho \in (0,1)$, the Markov process is said to be $\rho$-mixing conditioned on the invariant events.
````{prf:remark}
${\mathbb T}$ being a strong contraction on ${\mathcal N}$ limits intertemporal dependence of the Markov process $\{X_t\}$.
````

Let ${\mathbb I}$ be the identity operator. When the conditional expectation operator ${\mathbb T}$ is a strong contraction, the operator $({\mathbb I} - {\mathbb T})^{-1}$ is well defined, bounded on ${\mathcal N}$, and equal to the geometric sum:[^geometricSeries]

```{math}
\left({\mathbb I} - {\mathbb T}\right)^{-1} f(x) = \sum_{j=0}^\infty {\mathbb T}^j f(x) = \sum_{j=0}^\infty E \left[ f(X_{t+j}) \vert X_t = x \right].
```

[^geometricSeries]: The geometric series after the first equality sign is well defined under the weaker restriction that ${\mathbb T}^m$ is a strong contraction for some integer $m\geq 1$.
````{prf:example}
Consider the Markov chain setting of subsection [](subsec:chain) with a transition matrix ${\mathbb P}$.
A stationary density ${\bf q}$ is a nonnegative vector that satisfies
```{math}
{\bf q}' {\mathbb P} = {\bf q}'
```
and ${\bf q} \cdot \textbf{1}_n = 1 $.
If the only column eigenvector of ${\mathbb T}$ associated with a unit eigenvalue is constant over states $i$ for which ${\sf q}_i > 0$, then the process is ergodic.
If in addition, the only eigenvector of ${\mathbb P}$ that is associated with an eigenvalue that has a unit norm 
(the unit eigenvalue might be complex) is constant over states $i$ for which ${\sf q}_i > 0$, then ${\mathbb T}^m$ is a strong contraction for
some integer $m \geq 1$.[^gelfand]
This implies that the process is ergodic. It also rules out the presence of periodic components that can be forecast perfectly.

[^gelfand]: This follows from Gelfand's Theorem, which asserts the following. Let ${\mathcal N}$ be the $n-1$ dimensional
space of vectors that are orthogonal to ${\bf q}$. ${\mathbb T}$ maps ${\mathcal N}$ into itself.
The spectral radius of ${\mathbb T}$ restricted to ${\mathcal N}$ is the maximum of the 
absolute values of the eigenvalues.
Gelfand's Theorem asserts that the spectral radius governs the behavior as $m$ gets large of the decay factor of the ${\mathbb T}$ transformation applied $m$ times. Provided that the spectral radius is less than one,
 the strong contraction property prevails for any $\rho < 1$ that is larger than the spectral radius.
````

(sec:limitapprox)=
## Limits of Multi-Period Forecasts

When a Markov process is aperiodic, there are interesting situations in which
```{math}
:label: limitexp
\lim_{j \rightarrow \infty} {\mathbb T}^j f(x) =  {\sf r}
```
for some ${\sf r} \in {\mathbb R}$, where convergence is either pointwise in $x$ or in the ${\mathcal L}^2$ norm.
Limit {eq}`limitexp` asserts that long-run forecasts do not depend on the current Markov state. {cite:t}`meyntweedie` provide a comprehensive treatment of such convergence.
Let $Q$ be a stationary distribution. Then it is necessarily true that
```{math}
\int {\mathbb T}^j f(x) Q(dx)  = \int f(x) Q(dx)
```
for all $j$. Thus,
```{math}
{\sf r} = \int f(x) Q (dx),
```
so that the limiting forecast is necessarily the mathematical expectation of $f(x)$ under a stationary distribution.
Here we have assumed that the limit point is a number and not a random variable; we have not assumed that the stationary distribution is unique.

Notice that if {eq}`limitexp` is satisfied, then any function $f$ that satisfies
```{math}
{\mathbb T} f = f
```
is necessarily constant with probability one. Also, if $\int f(x) Q(dx) = 0$ and convergence is sufficiently fast, then
```{math}
:label: eq:sum_limit
\lim_{N \rightarrow \infty} \sum_{j=0}^N {\mathbb T}^j f(x)
```
is a well-defined function of the Markov state. We shall construct the limit in {eq}`eq:sum_limit` when we extract martingales from additive functionals in chapter [](chap:add). 




A set of sufficient conditions for the convergence outcome
```{math}
:label: converge
\lim_{j \rightarrow \infty} {\mathbb T}^j f (x^*) \rightarrow \int f(x) Q(dx)
```
for each $x^* \in {\mathcal X}$ and each bounded $f$ is:[^convergence_footnote]

(cond:suffice)= 
**Stability conditions:**


A Markov process with stationary distribution $Q$ satisfies:

(condi)= 
(i) For any $f \ge 0$ such that $\int f(x) Q(dx) > 0$, ${\mathbb M}_p f(x) > 0$ for all $x \in {\mathcal X}$ with $Q$ measure one and all positive integers $p \ge 0$, where the operator ${\mathbb M}_p$ is defined in equation not provided.

(fellercond)=
(ii) ${\mathbb T}$ maps bounded continuous functions into bounded continuous functions, i.e., the Markov process is said to satisfy the Feller property.

(condsupport)=
(iii) The support of $Q$ has a nonempty interior in ${\mathcal X}$.

(drift)=
(iv) ${\mathbb T} V(x) - V(x) \le -1$ outside a compact subset of ${\mathcal X}$ for some nonnegative function $V$.


[^convergence_footnote]: Restriction {eq}`converge` is stronger than ergodicity. It rules out periodic processes, although we know that periodic processes can be ergodic.

We encountered condition [(i)](condi) in our section [](sec:MarkErgodic). Condition [(iv)](drift) is a **drift condition** for stability that requires that we find a function $V$ that satisfies the requisite inequality. Heuristically, the drift condition says that outside a compact subset of the state space, application of the conditional expectation operator pushes the function inward. The choice of $-1$ as a comparison point is made only for convenience, since we can always multiply the function $V$ by a number greater than one. Thus, $-1$ could be replaced by any strictly negative number. In section [](sec:VAR44), we will apply conditions [(i)](condi) - [(iv)](drift)  to verify ergodicity of a vector autoregression.

(sec:VAR44)=
## Vector Autoregressions

A square matrix $\mathbb{A}$ is said to be *stable* when all of its eigenvalues have absolute values that are strictly less than one. For a stable $\mathbb{A}$, suppose that

```{math}
X_{t+1} = \mathbb{A} X_t + \mathbb{B} W_{t+1},
```
where $\{ W_{t+1} : t = 1,2,... \}$ is an i.i.d. sequence of multivariate normally distributed random vectors with mean vector zero and covariance matrix $I$ and that $X_0 \sim \mathcal{N}(\mu_0, \Sigma_0)$. This specification constitutes a first-order *vector autoregression*.

Let $\mu_t = E X_t$. Notice that

```{math}
\mu_{t+1} = \mathbb{A} \mu_t.
```
The mean $\mu$ of a stationary distribution satisfies

```{math}
:label: stationary_mean
\mu = \mathbb{A} \mu.
```
Because we have assumed that $\mathbb{A}$ is a stable matrix, $\mu =0$ is the only solution of {eq}`stationary_mean`, so the mean of the stationary distribution is $\mu = 0$.

Let $\Sigma_{t} = E(X_t - \mu_t) (X_t - \mu_t)'$ be the covariance matrix of $X_t$. Then

```{math}
\Sigma_{t+1} = \mathbb{A} \Sigma_t \mathbb{A}' + \mathbb{B}\mathbb{B}'.
```

For $\Sigma_t = \Sigma$ to be invariant over time, it must satisfy the discrete Lyapunov equation

```{math}
:label: eq:Sylvester
\Sigma = \mathbb{A} \Sigma \mathbb{A}' + \mathbb{B}\mathbb{B}'.
```

When $\mathbb{A}$ is a stable matrix, this equation has a unique solution for a positive semidefinite matrix $\Sigma$.

Suppose that $\Sigma_0 = 0$ (a matrix of zeros) and for $t \geq 1$ define the matrix

```{math}
\Sigma_t = \sum_{j=0}^{t-1} \mathbb{A}^j \mathbb{B}\mathbb{B}'(\mathbb{A}^j)'.
```
The limit of the sequence $\{\Sigma_t\}_{t=0}^{\infty}$ is

```{math}
\Sigma = \sum_{j=0}^{\infty} \mathbb{A}^j \mathbb{B}\mathbb{B}'(\mathbb{A}^j)',
```

which can be verified to satisfy Lyapunov equation {eq}`eq:Sylvester`. Thus, $\Sigma$ equals the covariance matrix of the stationary distribution.[^covmatrix] Similarly, for all $\mu_0 = E X_0$
```{math}
\mu_t = \mathbb{A}^t \mu_0,
```
converges to zero, the mean of the stationary distribution. The linear structure implies that the stationary distribution is Gaussian with mean $\mu$ and covariance matrix $\Sigma$. 

[^covmatrix]: To verify the asserted equality, notice that $\sum_{j=0}^\infty \mathbb{A}^j \mathbb{B}\mathbb{B}' \mathbb{A}^{j \prime} = \mathbb{A} ( \sum_{j=0}^\infty \mathbb{A}^j \mathbb{B}\mathbb{B}' \mathbb{A}^{j \prime} )\mathbb{A}' + \mathbb{B} \mathbb{B}'$. 

To verify ergodicity, we suppose that the covariance matrix $\Sigma$ of the stationary distribution has full rank and verify [Stability conditions](cond:suffice). Condition [(iii)](condsupport) since the covariance matrix has fully rank. Furthermore, $\Sigma_t$ has full rank for some $t$, which guarantees that the process is irreducible and aperiodic so that condition [(i)](condi) is satisfied. As a candidate for $V(x)$ in condition [(iv)](drift), take $V(x) = |x|^2$. Then

```{math}
\mathbb{T} V(x) = x'\mathbb{A}'\mathbb{A} x + \text{trace}(\mathbb{B}'\mathbb{B})
```

so

```{math}
\mathbb{T} V(x) - V(x) = x'(\mathbb{A}'\mathbb{A} - \mathbb{I})x + \text{trace}(\mathbb{B}'\mathbb{B}).
```

That $\mathbb{A}$ is a stable matrix implies that $\mathbb{A}'\mathbb{A} - \mathbb{I}$ is negative definite, so that drift condition [iv](drift)  is satisfied for $|x|$ sufficiently large.
The Feller property [(ii)](fellercond)  can also be verified by applying the The Dominated Convergence Theorem, since the functions used for the verification are bounded and continuous.  Thus, having checked [Stability conditions](cond:suffice), we have verified the ergodicity of the VAR.

## Estimating Vector Autoregressions

Let $Y_{t+1}$ be one of the entries of $X_{t+1}$, and consider the regression equation:

```{math}
Y_{t+1} = \beta \cdot X_{t} + U_{t+1},
```

where $U_{t+1}$ is a least squares residual. To express uncertainty about $\beta$ in the spirit of chapter [](chap:process), we allow it to be random. Letting ${\mathfrak J}$ be the set of invariant events, we presume that the *random vector* $\beta$ is measurable with respect to ${\mathfrak J}$, meaning that it is revealed by events in ${\mathfrak J}$.  For the case in which we have $n$ possibly different models, the invariant events reveal which of the models generates the data.

The first-order condition for minimizing the expected value of $U_{t+1}^2$ requires that the regression residual $U_{t+1}$ be orthogonal to $X_t$:

```{math}
E\left( X_t U_{t+1} \vert {\mathfrak J} \right) = 0.
```

Then

```{math}
:label: eq:LSorth101
E\left( X_{t}Y_{t+1} \vert {\mathfrak J} \right) = E\left[ X_{t} (X_{t})' \vert {\mathfrak J} \right] \beta ,
```

which uniquely pins down the regression coefficient $\beta$ provided that the matrix $E\left[ X_t (X_t)' \vert {\mathfrak J} \right]$ is nonsingular with probability one. Notice that

```{math}
{\frac 1 N}\sum_{t=1}^N  X_{t} Y_{t+1} \rightarrow E\left( X_{t} Y_{t+1} \vert {\mathfrak J} \right)  
```

```{math}
{\frac 1 N}\sum_{t=1}^N  X_{t} (X_{t})' \rightarrow E\left( X_{t} (X_{t})' \vert {\mathfrak J} \right),
```

where convergence is with probability one. Thus, from equation {eq}`eq:LSorth101` it follows that a consistent estimator of $\beta$ is a $b_N$ that satisfies

```{math}
{\frac 1 N}\sum_{t=1}^N  X_{t} Y_{t+1} = {\frac 1 N}\sum_{t=1}^N  X_{t} (X_{t})' b_N.
```

Solving for $b_N$ gives the familiar least squares formula:

```{math}
b_N = \left[\sum_{t=1}^N  X_{t} (X_{t})' \right]^{-1} \sum_{t=1}^N  X_{t} Y_{t+1}.
```

Note how statements about the consistency of $b_N$ are conditioned on ${\mathfrak J}$.  This conditioning is necessary when we do not know which among a family vector autoregressions generates the data.  

(sec:VAR_inf_past)=
## Inventing a Past Again

In section [](sec:inventing_past), we invented an infinite past for a stochastic process. Here we invent an infinite past for a vector autoregression in a way that is equivalent to drawing an initial condition $X_0$ at time $t=0$ from the stationary distribution ${\mathcal N}(0, \Sigma_\infty)$, where $\Sigma_\infty$ solves the discrete Lyapunov equation {eq}`eq:Sylvester`, namely, $\Sigma_\infty = {\mathbb A} \Sigma_\infty {\mathbb A}' + {\mathbb B} {\mathbb B}' $.

Thus, consider the vector autoregression
```{math}
X_{t+1} = {\mathbb A} X_t + {\mathbb B} W_{t+1}
```
where ${\mathbb A}$ is a stable matrix, $\{W_{t+1}\}_{t=-\infty}^\infty$ is now a two-sided infinite sequence of i.i.d. ${\mathcal N}(0,I)$ random vectors, and $t$ is an integer. We can solve this difference equation backwards to get the moving average representation

```{math}
X_{t} = \sum_{j=0}^\infty {\mathbb A}^j {\mathbb B} W_{t -j} .
```
Then
```{math}
E\left[X_t \left(X_t \right)'\right] = \sum_{j=0}^\infty {\mathbb A}^j {\mathbb B} {\mathbb B}' \left( {\mathbb A}^j \right)' = \Sigma_\infty
```
where  $\Sigma_\infty$ is also the unique positive semidefinite matrix that solves $\Sigma_\infty = {\mathbb A} \Sigma_\infty {\mathbb A}' + {\mathbb B} {\mathbb B}'$.
