(chap:book)=
# GMM Estimation [^note1]
[^note1]: These notes are very preliminary.

Related papers: 
1. [Hurwicz (1966)](https://www.sciencedirect.com/science/article/abs/pii/S0049237X09705907): On the Structural Form of Interdependent Systems

---


## Introduction

Generalized method moments (GMM) estimation studies a family of estimators constructed from partially specified or partially misspecified models. Direct application of likelihood methods sometimes can be challenging to construct, and  GMM methods may be tractable alternatives. By studying an entire family of estimators, we are able to make relative accuracy comparisons among the entire family.

This paper takes statistical consistency as given. Supporting arguments for this chapter can be obtained with direct extensions of the Law of Large Numbers. Such extensions often entail Laws of Large numbers applied to so-called random processes (indexed by say a parameter vector of interest) instead of a random vector.

Throughout this chapter, we will condition on invariant events even though we will suppress this dependence when we write conditional expectations.  Given the partially specified or misspecified nature of the model, much more than a simple parameter vector is reflected by this conditioning.   


## Formulation

We study a family of GMM estimators of an unknown parameter vector $\beta$ constructed from theoretical restrictions on conditional or unconditional moments of functions $\phi$ that depend on $\beta$ and on a random vector $X_t$ that is observable to an econometrician.

As a starting point, we consider a class of restrictions large enough to include examples of both conditional and unconditional moment restrictions. Members of this class take the form 

```{math}
:label: GMMpopulation
E \left[ {A_t}' \phi(X_t, b) \right] = 0  \textrm{ if and only if } b = \beta
``` 

for all sequences of selection matrices $A  \in \mathcal A$ where $A = \{A_t : t \ge 1\} $ and where

- the vector of functions $\phi$ is $r$ dimensional.
- the unknown parameter vector $\beta$ is $k$ dimensional.
- ${\mathcal A}$ is a collection of sequences of (possibly random) selection matrices that characterize valid moment restrictions.
- $A_t$ denotes a time $t$ selection matrix for a subset of the valid moment restrictions that is used to construct a particular statistical estimator $b$ of $\beta$.
- the mathematical expectation is taken with respect to a statistical model that generates the $\{X_t : t \ge 1 \}$ process (captured implicitly by conditioning on invariant events).  

A sample counterpart of the population moment conditions {eq}`GMMpopulation` is 

```{math}
:label: GMMsample
\frac 1 N \sum_{t=1}^N {A_t}'  \phi(X_t, b_N) = 0 .
```
Applying a Law of Large Numbers to {eq}`GMMsample` motivates a "generalized method of moments" estimator $b_N$ of the $k \times 1$ vector $\beta$.

Different sequences of selection matrices $\{A_t : t \ge 1 \}$ and $\{\widetilde A_t: t \ge 1\}$ generally give rise to different properties for the estimator $b_N$.  An exception is when 

```{math}
{\widetilde A}_t =  A_t{\mathbb L} 
```

for some $k \times k$ nonsingular matrix ${\mathbb L}$ do distinct selection matrices ${\widetilde A}_t$ and $A_t$ give rise to the same $b_N$.

We study limiting properties of estimator $b_N$ conditioned on a statistical model. In many settings, the parameter vector $\beta$ only incompletely characterizes the statistical model. In such settings, we are led in effect to implement a version of what is known as *semi-parametric* estimation: while $\beta$ is the finite-dimensional parameter vector that we want to estimate, we acknowledge that, in addition to $\beta$, a potentially infinite dimensional nuisance parameter vector pins down the complete statistical model on which we condition when we apply the law of large numbers and other limit theorems.

````{prf:example}
:label: ex:unconditional

**Unconditional moment restrictions**  

Suppose that 

```{math}
E \left[ \phi(X_t, \beta) \right] = 0
```

where $r \ge k$. Let $\mathcal{A}_t$ be the set of all constant $r \times k$ matrices $\mathbb{A}$ of constants. Rewrite the restrictions as:

```{math}
{\mathbb{A}}' E \left[ \phi(X_t, \beta) \right] = 0
```

for all $r \times k$ matrices $\mathbb{A}$. {cite}`sargan58` and {cite}`hansen82` assumed moment restrictions like these.  For instance, 

```{math}
\phi(X_t,\beta) = Z_t \eta(Y_t, \beta )'
```
where $Z_t$ is an $r$ dimensional vector of *instrumental variables* and $\eta(Y_t, \beta)$ is a scalar disturbance term in an equation of interest.  The vector of *instrumental variables* are presumed to be uncorrelated with $Z_t$. 

````

````{prf:example}
:label: ex:conditional
**Conditional moment restrictions**

Assume the conditional moment restrictions 
```{math}
E\left[\phi(X_{t}, \beta) \mid {\mathfrak A}_{t-\ell} \right] = 0
```
for a particular $\ell \ge 1$ and $Y_t = X_t$. Let $\mathcal A_{t}$ be the set of all $r \times k$ matrices, $A_t$, of bounded random variables that are ${\mathfrak A}_{t-\ell}$ measurable.  
Then the preceding conditional moment restrictions are mathematically equivalent to the unconditional moment restrictions
```{math}
E\left[{A_t}'\phi(Y_{t}, \beta) \right] = 0
```
for all random matrices $A_t \in {\mathcal A}_t$. This formulation is due to {cite}`hansen85`. Also see a closely analysis of {cite}`chamberlain87`.
````

Collections $\mathcal A$ of selection processes for both of these examples satisfy the following "linearity" restriction.


(res1)=
**Restriction 9.1.** If $A^1$ and $A^2$ are both in $\mathcal A$ and $\mathbb{L}_1$ and $\mathbb{L}_2$ are $k \times k$ matrices of real numbers, then $A^1 \mathbb{L}_1 + A^2\mathbb{L}_2$ is in $\mathcal A$.

A common practice is to use the idea provided in {prf:ref}`ex:conditional` while substantially restricting the set of moment conditions used for parameter estimation. Thus, from a collection of conditional moment restrictions, we can create unconditional moment restrictions like those in {prf:ref}`ex:unconditional` and thereby reduce the class of GMM estimators under consideration. For instance, let $A_t^1$ and $A_t^2$ be two *ad hoc* choices of selection matrices. Form
```{math}
\phi^+(X_t, b) = \begin{bmatrix} {A_t^1}' \\ {A_t^2}'  \end{bmatrix} \phi(X_t, b)
```
where $X_t$ now includes variables used to construct $A_t^1$ and $A_t^2$. We presume that no linear combination of columns of $A_t^2$ duplicate any columns in $A_t^1$. Otherwise, we would omit such columns and adjust $\phi^+$ accordingly. Let $r^+ \ge r$ denote the remaining non-redundant columns. We use $r^+ \times k$ selection matrices $\mathbb A$ to form moment conditions 
```{math}
{\mathbb A}'E\left( \phi^+(X_t, b)\right) = 0
```
and study an associated family of GMM estimators. This strategy reduces an infinite number of moment conditions to a finite number. There are extensions of this approach. For instance, we could use more than two $A_t^j$'s to construct $\phi^+$.

````{prf:remark}
:label: momentmatch1

"Moment matching" estimators are another special case of {prf:ref}`ex:unconditional`.
Suppose that

```{math}
\phi(X_t, b) = \psi(X_t) - \kappa(b)
```

where

```{math}
E\left[\psi(X_t) \right] = \kappa(\beta).
```

The random vector  $\psi(Y)$ defines moments to be matched and $\kappa(\beta)$ are population values of those moments under a statistical model
with parameter vector $\beta$. Often that statistical model is a “structural” economic model with nonlinearities and other complications that, for a given, $\beta$ make it challenging to compute the moments $E\left[\psi(X_t) \right]$ analytically. To proceed, the proposal is to approximate those moments for a given $b$ by computing a sample mean from a long simulation of the statistical model at parameter vector $b$. By running simulations and computing associated sample means for many alternative $b$ vectors, we can assemble an approximation to the function $\kappa(b)$. {cite}`LeeIngram:1991` and {cite}`DuffieSingleton:1993` used versions of this approach. Notice that in contrast to some other applications of GMM estimation that allow the appearance of unknown nuisance parameters in the statistical model assumed to generate the data, this approach assumes that, given $b$, the model completely determines a sample path that we can at least simulate.

````

````{prf:remark}
"Indirect inference" works with two statistical models in hand: (1) a "structural economic model" with a vector of parameters $\beta$ that characterize preferences, technology, information flows, and other features of the theoretical economic model; and (2) an "auxiliary model" with no pretense of being "structural" in terms of economic theory and having a vector of parameters $\delta$ that let the model fit well. Although the structural model can be solved and simulated on a computer, it is too complicated to allow writing down its likelihood process analytically. The likelihood process for the auxiliary model can be calculated analytically. "Moments matching" estimation in the style of {cite}`GallantTauchen:1996` proceeds in two steps, the first of which uses maximum likelihood estimation of parameters of the auxiliary model prepares a random vector $\psi(X_t)$ whose moments are to be matched, the second of which proceeds as in {prf:ref}`momentmatch1` to use simulations of the structural model to approximate the function $\kappa (\beta)$. In the first step, the parameter vector $\delta$ of the auxiliary model is estimated by the method of maximum likelihood and the sample path of the associated score vector is evaluated at the maximum likelihood estimate $\hat{\delta}$. As an input into the second step, the associated Fisher information matrix is computed. The second step forms a GMM criterion consisting of a quadratic form in the score vector with weighting matrix being the inverse of the Fisher information matrix computed in the first step. Repeated simulations of the structural model are used to search for a $b_N$ that best matches score-criterion from the auxiliary model.

[^orphan-footnote-with-references]: Important related approaches use a misspecified maximum likelihood ({cite}`Smith:1993` and {cite}`GourierouxMonfortRenault:1993`) or the score increment of such a likelihood ({cite}`GallantTauchen:1996`) to summarize empirical evidence and use model simulation to account for the misspecification. The method is applied either for reasons of computational simplicity or because the researcher wants to feature believed to be robust to misspecification of the statistical model. _Lars -- I'd like to discuss how the ideas in this footnote are phrased. "misspecified" is a loaded term here (and in the text above) -- susceptible to multiple interpretations in the ears of a reader I fear._
````
## Central limit approximation

The process 
```{math}
\left\{ \sum_{t=1}^N{A_t}'  \phi(X_t, \beta) : N \ge 1\right\} .
```
can be verified to have stationary and ergodic increments conditioned on the statistical model. So there exists a Proposition 2.2.2 decomposition of the process. Provided that
```{math}
E\left[  {A_t}'  \phi(X_t, \beta)\right] =0 
```
under the statistical model that generates the data, the trend term in the decomposition of Proposition 2.2.2 is zero, implying that the martingale dominates the behavior of sample averages for large $N$. In particular, Proposition 2.3.1 in Chapter 2 gives a central limit approximation for
```{math}
\frac 1 {\sqrt{N}}  \sum_{t=1}^N{A_t}'  \phi(X_t, \beta) 
```
Let $A = \{ A_t : t \ge 0 \}$ and suppose that
```{math}
E \left[ \sum_{j=0}^\infty {A_{t+j}}' \phi(X_{t+j}, \beta) \mid {\mathfrak A}_t \right] 
```
converges in mean square. Define the one-step-ahead forecast error:
```{math}
G_t(A) = E \left[ \sum_{j=0}^\infty {A_{t+j}}' \phi(X_{t+j}, \beta) \mid {\mathfrak A}_t \right] - E \left[ \sum_{j=0}^\infty {A_{t+j}}' \phi(X_{t+j}, \beta) \mid {\mathfrak A}_{t-1} \right] 
```
Paralleling the construction of the martingale increment in Proposition 2.2.2,
```{math}
\frac 1 {\sqrt{N}}  \sum_{t=1}^N{A_t}'  \phi(X_t, \beta) \approx {\frac 1 {\sqrt{N} }} \sum_{t=1}^N G_t(A) 
```
where by the approximation sign $\approx$ we intend to assert that the difference between the right side and left side converges in mean square to zero as $N \rightarrow \infty$. Consequently, the covariance matrix in the central limit approximation is
```{math}
E \left[ G_t(A) G_t(A)' \right].
```

Recall [Restriction 9.1](res1). For the preceding construction of the martingale increment, it is straightforward to verify that
```{math}
G_t( A^1 {\mathbb L}_1 +  A^2{\mathbb L}_2) = ({\mathbb L}_1)' G_t(A^1) + ({\mathbb L}_1)'G_t(A^2)
```
follows from the linearity of conditional expectations.

````{prf:example}
Consider again {prf:ref}`ex:unconditional` in which  $A_t = \mathbb{A}$ for all $t\ge 0$ and 
```{math}
G_t(A) = \mathbb{A}' F_t
```
where
```{math}
F_t =  
E \left[ \sum_{j=0}^\infty  \phi(X_{t+j}, \beta) \mid \mathfrak{A}_t \right] - E \left[ \sum_{j=0}^\infty \phi(X_{t+j}, \beta) \mid \mathfrak{A}_{t-1} \right]. 
```
Define the covariance matrix
```{math}
\mathbb{V} = E \left( F_t F_t ' \right)
```
and note that 
```{math}
E\left[ G_t(A)G_t(A)' \right] = \mathbb{A}' \mathbb{V} \mathbb{A}  .
```
````

````{prf:example}
In {prf:ref}`ex:conditional`

```{math}
E\left[\phi(Y_{t}, \beta) \mid {\mathfrak A}_{t-\ell} \right] = 0
```

and hence

```{math}
E\left[{A_t}' \phi(Y_{t}, \beta) \mid {\mathfrak A}_{t-\ell} \right] = 0
```

whenever entries of $A_t$ are restricted to be ${\mathfrak A}_{t-\ell}$ measurable.  
Consequently

```{math}
E\left[{A_{t+j}}' \phi(Y_{t+j}, \beta) \mid {\mathfrak A}_{t} \right] = 0
```

for $j \ge \ell$ so that the infinite sums used to construct $G_t(A)$ simplify to finite sums.
````
## Mean value approximation

Write

```{math}
\frac{1}{\sqrt{N}}
\sum_{t=1}^NA_t' \phi(X_t, b_N)  & \approx
\frac{1}{\sqrt{N}} \sum_{t=1}^NA_t' \phi(X_t, \beta) +
\frac{1}{N} \sum_{t=1}^NA_t' \left[\frac{\partial \phi}{\partial b'} (X_t, \beta)\right] \sqrt{N} (b_N - \beta) \\
& \approx  \frac{1}{\sqrt{N}} \sum_{t=1}^NA_t' \phi(X_t, \beta)  + \nabla(A)'\sqrt{N} (b_N - \beta) 
```

where

```{math}
\nabla(A) \overset{\text{def}}{=} E\left(\left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right]' A_t \right).
```

Since 
```{math}
\frac{1}{\sqrt{N}}
\sum_{t=1}^NA_t' \phi(X_t, b_N) \approx 0, 
```
```{math}
\nabla(A)' \sqrt{N} (b_N - \beta)  \approx  - \frac{1}{\sqrt{N}} \sum_{t=1}^NA_t' \phi(X_t, \beta) .
```

So long as $\nabla(A)$ is nonsingular,

```{math}
\sqrt{N} (b_N - \beta) \approx - \left[\nabla(A)'\right]^{-1}  \frac{1}{\sqrt{N}} \sum_{t=1}^NA_t' \phi(X_t, \beta) .
```

This approximation underlies an "efficiency bound" for GMM estimation. Notice that the covariance matrix in a central limit approximation is:

```{math}
\textbf{cov}(A) = \left[\nabla(A)'\right]^{-1} E\left[ G_t(A) G_t(A)' \right] \left[\nabla(A) \right]^{-1}
```

We want to know how small we can make this matrix by choosing a selection process.

````{prf:example}
Consider again {prf:ref}`ex:unconditional`.  In this case $A_t = {\mathbb A}$ for all $t \ge 0$ and 
```{math}
\nabla(A) =  {\mathbb D}' {\mathbb A}  
```
where 
```{math}
{\mathbb D} \overset{\text{def}}{=} E\left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right]  
```
and
```{math}
\textrm{ cov}({\mathbb A} ) = \left({\mathbb A}' {\mathbb D}\right)^{-1}  {\mathbb A}' {\mathbb V}  {\mathbb A}  \left({\mathbb D}' {\mathbb A}\right)^{-1}
```
````

For purposes of devising a test of the "over-identifying restrictions," let $B = \{ B_t : t \ge 0\}$ be an $r  \times {\tilde k}$ matrix process constructed to verify
```{math}
E  \left[{B_t}'  \phi(X_t, \beta) \right] = 0.
```


(res2)=
**Restriction 9.2.** For any ${\tilde k} \times k$ matrix of real numbers ${\mathbb K}$, $ B {\mathbb K} \in {\mathcal A}$.   


Thus, we can build selection processes for estimation equations from the columns of the process $B$.  

Suppose that
```{math}
E \left[ \sum_{j=0}^\infty {B_{t+j}}' \phi(X_{t+j}, \beta) \mid {\mathfrak A}_t \right] 
```
converges in mean square so that we can apply a central limit approximation.  
Construct
```{math}
\widetilde {\nabla}(B)  \overset{\text{def}}{=} E\left( \left[\frac {\partial \phi}{\partial b'} 
(X_t, \beta)\right]' B_t \right) .
```
Since [Restriction 9.2](res2) is satisfied, notice that
```{math}
\widetilde {\nabla}(B) {\mathbb K} = \nabla (B {\mathbb K} )
```
for all ${\tilde k} \times k$ matrices ${\mathbb K}$ of real numbers.  
 
By imitating an earlier argument
```{math}
\begin{align*}
\frac 1 {\sqrt{N}} 
\sum_{t=1}^N {B_t}' \phi(X_t, b_N) \approx &
  \frac 1 {\sqrt N} \sum_{t=1}^N {B_t}' \phi(X_t, \beta)  + \widetilde {\nabla}(B)' \sqrt{N} (b_N - \beta) \cr
  \approx &   \frac 1 {\sqrt N} \sum_{t=1}^N {B_t}' \phi(X_t, \beta)\cr &  - \widetilde {\nabla}(B)' \nabla(A)^{-1}  \frac 1 {\sqrt N} \sum_{t=1}^N {A_t}' \phi(X_t, \beta) \cr
\approx & \frac 1 {\sqrt N} \sum_{t=1}^N \left[{B_t}'  - \widetilde {\nabla}(B) '\left[\nabla(A)'\right]^{-1} {A_t}'\right] \phi(X_t, \beta)  
\end{align*} 
```
Notice that if $A_t = B_t$, then the right side is zero and the limiting distribution is degenerate.  This approximation is used to construct tests that account for having used GMM to estimate a parameter vector $\beta$.

````{prf:example}
Consider again unconditional moment restrictions specified in {prf:ref}`ex:unconditional`. Let the selection process for testing be constant over time so that
$B_t = {\mathbb B}$. Then 

```{math}
\frac 1 {\sqrt{N}} 
\sum_{t=1}^N {B_t}' \phi(X_t, b_N) \approx  \frac 1 {\sqrt N} \sum_{t=1}^N 
\left[{\mathbb B}'  - {\mathbb B}' {\mathbb D} \left({\mathbb A}' {\mathbb D} \right)^{-1} {\mathbb A}'\right]\phi(X_t, \beta) .
```
````
## GMM Efficiency Bound

Recall
```{math}
\textrm{cov}(A) = \left[\nabla(A)'\right]^{-1} E\left[ G_t(A) G_t(A)' \right] \left[\nabla (A) \right]^{-1}
```
We seek a greatest lower bound on the covariance matrix on the right.

1. Suppose that $\left[\nabla(A)'\right]^{-1}$ is nonsingular and impose that 
   ```{math}
   \left[\nabla(A)\right] = {\mathbb I}
   ```
   If not, post multiply $A$ by a nonsingular matrix ${\mathbb J}$. That leaves the GMM estimator unaltered.
   Thus, we have 
   ```{math}
   \textrm{ cov}(A) =  E\left[ G_t(A) G_t(A)' \right] 
   ```
   subject to $\left[\nabla(A)\right] = {\mathbb I}$

2. Find an $A^d$ such that for all $A \in {\mathcal A}$
   ```{math}
   :label: GMM_first
   \nabla(A) = E\left[ G_t(A^d) G_t(A)' \right]  .
   ```

3. Form
   ```{math}
   A_t^* = A^d_t \left( E\left[ G_t(A^d) G_t(A^d)' \right]\right)^{-1} 
   ```
   for all $A \in \mathcal A$. These form a set of first-order sufficient conditions for our constrained  minimization problem.  
   Then
   ```{math}
   G_t(A^*) = \left( E\left[ G_t(A^d) G_t(A^d)' \right]\right)^{-1} G_t(A^d)
   ```
   and
   ```{math}
   E \left[ G_t(A^*) G_t(A)'\right] = \left( E\left[ G_t(A^d) G_t(A^d)' \right]\right)^{-1}
   ```
   provided that $\left[\nabla(A)\right] = {\mathbb I}.$  

4. Therefore,
   ```{math}
   0 \leq E \left( \left[ G_t(A) - G_t(A^*) \right] \left[ G_t(A) - G_t(A^*) \right]' \right) = \textrm{cov}(A) - \left( E\left[ G_t(A^d) G_t(A^d)' \right]\right)^{-1} .
   ```



(res:eff_bound)=
**Result 9.1.** Given a solution to equation {eq}`GMM_first`
```{math}
:label: eq:effbound105
\inf_{A \in \mathcal A}   \textrm{cov}(A) = \left( E\left[ G_t(A^d) G_t(A^d)' \right]\right)^{-1}
```


````{prf:remark}
In the result [9.1](res:eff_bound) efficiency bound, we might be tempted to think that $G_t(A^d)$ plays the same role that the "score vector" increment does in maximum likelihood estimation. But because there is a possibly infinite dimensional vector of nuisance parameters here, a better analogy is that $G_t(A^d)$ acts much like the residual vector in a regression of parameters of interest score increments on nuisance parameter score increments. 
By undertaking to infer the parameter vector $\beta$ from conditional or unconditional moment restrictions, we have purposefully pushed all nuisance parameters into the background.
````

````{prf:remark}
The representation  

```{math}
E\left(\left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right]' A_t \right) = \nabla(A) = E\left[ G_t(A^d) G_t(A)' \right]  
```

used to compute the efficiency bound is an application of the Riesz Representation Theorem. To understand this, introduce the $k$-dimensional coordinate vectors ${\sf u}_i$ for $i=1,2,...,k$ and consider:

```{math}
:label: Rrep
({\sf u}_i)'\nabla(A) {\sf u}_j = E\left(\left[\frac {\partial \phi}{\partial b_i} (X_t, \beta)\right] \cdot  A_t {\sf u}_j\right) .
```

Note that

- The integer $i$ selects the coordinate of $b$ with respect to which we are differentiating.
  
- If $A^1$ and $A^2$ are both in ${\mathcal A}$, then so are linear combinations. Therefore $({\sf u}_i)'\nabla(A) {\sf u}_j$ is a linear functional defined on a linear space of random variables of the form $\phi(X_t, \beta)' A_t {\sf u}_j$ for a given $i$.

- The martingale approximations for the scalar process with stationary increments 

  ```{math}
  \left\{ \sum_{t=1}^N  \phi(X_t, \beta)' A_t {\sf u}_j : N \ge 1 \right\}
  ```
  
  has martingale increment $G_t(A) {\sf u}_j$.

- The Riesz Representation Theorem asserts that the linear functional $({\sf u}_i)'\nabla(A) {\sf u}_j$ can be represented as an inner product 

  ```{math}
  ({\sf u}_i)'\nabla(A) {\sf u}_j = E\left[ R_t G_t(A) {\sf u}_j \right]
  ```
  
  where the scalar random variable $R_t$ is in the mean square closure of

  ```{math}
  \left\{ G_t(A) {\sf u}_j : A \in {\mathcal A} \right\}.
  ```

- We can represent $R_t$ as

  ```{math}
  R_t = G_t(A^d){\sf u}_j
  ```
  
  for some selection process $A_d \in {\mathcal A}$ or more generally as a limit point of a sequence of such selection processes.

The preceding construction pins down row $j$ of $A^d$. Repeating an analogous construction for each $j = 1,2,...,k$ gives the selection matrix $A^d$.

The GMM efficiency bound presumed that we could solve equation {eq}`Rrep`. The Riesz Representation Theorem requires that $R_t$ be in a mean square closure of a linear space. Provided that the linear functionals $({\sf u}_i)'\nabla(A) {\sf u}_j$ are mean square continuous, the efficiency bound can be represented in terms of the limit point of a sequence of GMM estimators associated with alternative selection processes even when the limit points are not attained.

````

````{prf:example}
Consider {prf:ref}`ex:unconditional` in which we assumed that $A_t = \mathbb{A}$.  Then

```{math}
\mathbb{A}' \mathbb{V} \mathbb{A}^d =  \mathbb{A}'\mathbb{D} .
```
Therefore,

```{math}
\mathbb{A}^d =  \mathbb{V}^{-1} \mathbb{D} 
```

and the GMM efficiency bound is

```{math}
\left(\mathbb{D}'  \mathbb{V}^{-1} \mathbb{D}\right)^{-1} .  
```
````

````{prf:example}
:label: ex:exampleeffb103
Consider again {prf:ref}`ex:conditional` in the special case in which $\ell = 1$.  
Let 
```{math}
E \left[ \phi(X_t, \beta) \phi(X_t, \beta)' \mid {\mathfrak A}_{t-1} \right] = V_{t-1} 
```
wish to solve the following equation for $A_t^d$
```{math}
:label: first_order
E\left( {A_t^d}' V_{t-1} A_t \right) = \nabla(A) = 
E\left(\left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right]' A_t \right) .
```
Given the flexibility in the choice of the random $A_t$ with entries that are ${\mathcal A}_{t-1}$ measurable, this equation is equivalent to
```{math}
V_{t-1} A_t^d = E\left( \left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right] \mid {\mathfrak A}_{t-1} \right)
```
where we have taken transposes of the expressions in {eq}`first_order`. Thus
```{math}
A_t^d = \left(V_{t-1}\right)^{-1} E\left( \left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right] \mid {\mathfrak A}_{t-1} \right) 
```
and the efficiency bound is:
```{math}
\left[  E\left( \left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right]' \mid {\mathfrak A}_{t-1} \right) \left(V_{t-1}\right)^{-1}E\left( \left[\frac {\partial \phi}{\partial b'} (X_t, \beta)\right] \mid {\mathfrak A}_{t-1} \right) \right]^{-1}.
```
````

````{prf:example}
:label: ex:2sls
**Two-stage least squares**. Add the following special restrictions to {prf:ref}`ex:exampleeffb103`. Suppose that $r=1$ and that 
$V_{t-1} = \mathsf{v} > 0$ where $\mathsf{v}$ is constant. Further suppose that

```{math}
\phi(X_t, b) = Y_t^1 - Y_t^2 \cdot b
```

Finally, suppose that

```{math}
E\left( Y_t^2 \mid \mathfrak{A}_{t-1} \right) = \Pi Z_{t-1}
```

where $Z_{t-1}$ has more entries than $Y_t^2$.
Notice that $\Pi$ can be computed as a least squares regression. Then

```{math}
A_t^d = \left(\frac{1}{\mathsf{v}}\right) {Z_{t-1}}' \Pi'
```

The scaling by $\frac{1}{\mathsf{v}}$ is inconsequential to the construction of a selection process. The matrix of regression coefficients can be replaced by the finite sample least squares regression coefficients without altering the statistical efficiency.
````

````{prf:remark}
{prf:ref}`ex:2sls` has a special structure that does not prevail in some important applications. For instance, suppose that $V_{t-1}$ depends on conditioning information so that a form of conditional heteroskedasticity is present. That dependence shows up in essential ways in how $A^d_t$ should be constructed. Further, suppose that the expectation $E\left( X_t^2 \mid {\mathfrak A}_{t-1} \right)$ potentially depends nonlinearly on $Z_{t-1}$. In that case, to attain or to approximate the efficiency bound, a least squares regression should account for potential nonlinearity. Finally, suppose that $\ell > 1$. Then even if the covariance structure is homoskedastic and conditional expectations are linear, the two-squares least square approach will no longer be statistically efficient. We again have to deploy
an appropriate martingale central limit approximation. In these circumstances, simply by mapping into the framework of {prf:ref}`ex:unconditional`, we can improve efficiency relative to least squares or two-stage least squares, for instance, by letting
```{math}
\phi(X_t, b) =  Z_{t-\ell}\left[Y_t^1 - \left(Y_t^2 \right)' b \right]
```
{cite}`HansenSingleton:1996` construct the efficiency bound in {prf:ref}`ex:conditional` for a linear data generating process.
````
## Statistical tests

First suppose that we have statistically efficient selection process. Recall the approximation

```{math}
\frac 1 {\sqrt{N}} 
\sum_{t=1}^N {B_t}' \phi(X_t, b_N) \approx \frac 1 {\sqrt N} \sum_{t=1}^N \left[{B_t}'  - \widetilde {\nabla}(B) '\left[\nabla(A^d)'\right]^{-1} {A_t^d}'\right] \phi(X_t, \beta)   .
```

Let ${\widetilde G}_t(B)$ denote the increment in the martingale approximation for

```{math}
\sum_{t=1}^N {B_t}' \phi(X_t, \beta) .
```

From the restrictions that we have imposed on the process $B$ used for constructing tests

```{math}
\widetilde {\nabla}(B) =  E\left[ G_t(A^d) G_t(B)' \right] .
```

Using both of these representations:

```{math}
:label: increment_regress
\frac 1 {\sqrt{N}} 
\sum_{t=1}^N {B_t}' \phi(X_t, b_N)  \approx \frac 1 {\sqrt N} \sum_{t=1}^N  {\widehat G}_t(B)
```

where

```{math}
{\widehat G}_t(B) \overset{\text{def}}{=} 
{\widetilde G}_t(B) - 
E\left[{\widetilde G}_t(B) G_t(A^d)'\right] \left( E\left[{ G}_t(A^d) G_t(A^d)'\right] \right)^{-1} G_t(A^d)     
```

The term, ${\widehat G}_t(B)$, that appears inside the sum on the right side of {eq}`increment_regress` is the population least squares residual from regressing ${\widetilde G}_t(B)$ onto $G_t(A^d)$. This regression residual can also be interpreted as a martingale increment for a stationary increments process.

Suppose that ${\widehat G}_t(B)$ has a nonsingular covariance matrix. Consider the quadratic form used for building a test:

```{math}
{\frac 1 N} \left[\sum_{t=1}^N \phi(X_t, b_N)' B_t \right] \left(E \left[{\widehat G}_t(B){\widehat G}_t(B)' \right]\right)^{-1} \left[\sum_{t=1}^N {B_t}' \phi(X_t, b_N)\right] \Rightarrow \chi^2 ({\tilde k} ) .
```

This test can be implemented in practice by replacing $E \left[{\widehat G}_t(B){\widehat G}_t(B)' \right]$ with a statistically consistent estimator of it. There is an equivalent way to represent this quadratic form:

```{math}
\begin{align}
{\frac 1 N} \sum_{t=1}^N \phi(X_t,b_N)'\begin{bmatrix} {B_t} &  {A_t^d} \end{bmatrix} & 
\left[E \left( \begin{bmatrix} {\widetilde G}_t(B) \cr G_t(A^d) \end{bmatrix} 
 \begin{bmatrix} {\widetilde G}_t(B)' &  G_t(A^d)' \end{bmatrix} \right)\right]^{-1} \\ 
 
& \left[\sum_{t=1}^N \begin{bmatrix} {B_t}' \cr   {A_t^d}' \end{bmatrix} \phi(X_t, b_N)\right] 
\end{align}
```

This equivalence follows because the inverse of the covariance matrix for the regression error ${\widehat G}_t(B)$ is the upper diagonal block of the inverse of the covariance matrix:

```{math}
E \left( \begin{bmatrix} {\widetilde G}_t(B) \cr G_t(A^d) \end{bmatrix} 
 \begin{bmatrix} {\widetilde G}_t(B)' &  G_t(A^d)' \end{bmatrix} \right)
```

````{prf:example}
Consider {prf:ref}`ex:unconditional` again. We have already shown that

```{math}
{\mathbb A}^d = {\mathbb V}^{-1} {\mathbb D} .
```

Suppose that we choose ${\mathbb B}$ with dimension $r \times (r-k)$ so that

```{math}
\begin{bmatrix} {\mathbb A}^d & {\mathbb B} \end{bmatrix}
```

has full rank. Then

```{math}
{\frac 1 N} \sum_{t=1}^N \phi(X_t,b_N)' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,b_N)' \Rightarrow \chi^2(r - k) .
```

If we replace $b_N$ with $\beta$ on the left side of the above limit we find

```{math}
{\frac 1 N} \sum_{t=1}^N \phi(X_t,\beta)' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,\beta)' \Rightarrow \chi^2(r)
```

The difference in the resulting $\chi^2$ distribution emerges because estimating $k$ free parameters reduces degrees of freedom by $k$. It is straightforward to show that

```{math}
{\frac 1 N} \sum_{t=1}^N \phi(X_t,\beta)' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,\beta)' - {\frac 1 N} \sum_{t=1}^N \phi(X_t,b_N)' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,b_N)' \Rightarrow \chi^2(k) ,
```

an approximation that is useful for constructing confidence sets for GMM estimates of parameter vector $\beta$.

````

````{prf:remark}
:label: rem:old5.2
To continue our study of {prf:ref}`ex:unconditional`, form the population problem:

```{math}
\min_b E\left[ \phi(X_t, b) \right]' {\mathbb V}^{-1}  E \left[\phi(X_t,b)\right] .
```
This has a minimizer at $b = \beta$ provided that the unconditional moment conditions are satisfied.  If $b = \beta$ is the only possible parameter vector that satisfies the population moment conditions, then $b = \beta$ is the unique solution to the population minimization problem stated here.  Suppose that we construct an estimator by solving a minimization problem:

```{math}
:label: GMM_min
\min_b  {\frac 1 N} \sum_{t=1}^N \phi(X_t,b)' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,b) .
```
First-order necessary conditions are 

```{math}
{\frac 1 N} \sum_{t=1}^N \left[\frac {\partial \phi }{\partial b'}(X_t,b_N)  \right]' {\mathbb V}^{-1} \sum_{t=1}^N \phi(X_t,b_N) = 0 .
```
Assume that we already know that  the solution $b_N$ of the above first-order conditions provides a consistent estimator of parameter vector $\beta$. Then we can show that 

```{math}
{\frac 1 N} \sum_{t=1}^N \left[\frac {\partial \phi }{\partial b'} (X,b_N)\right]  \rightarrow {\mathbb D}
```
where convergence is with probability one.  Thus, in this case  the selection matrix 

```{math}
{\mathbb A} = {\mathbb V}^{-1} {\mathbb D} 
```
provides an estimator that attains the efficiency bound.  The limiting distribution of the minimizer of criterion {eq}`GMM_min` is $\chi^2$ with $r - k$ degrees of freedom.
````

````{prf:remark}
There is an interesting variation of the approach described in {prf:ref}`rem:old5.2`. For any $b$, let $\mathbb{V}(b)$ be the population covariance matrix in the martingale increment used in the Central Limit approximation for the process

```{math}
\frac{1}{\sqrt{N}} \sum_{t=1}^N \phi(X_t, \beta) .
```

Assume that $\mathbb{V}(b)$ is nonsingular for every $b$ in a parameter space. Form the population minimization problem:

```{math}
\min_b E\left[ \phi(X_t, b) \right]' \left[\mathbb{V}(b)\right]^{-1} E \left[\phi(X_t,b)\right] .
```

If $b = \beta$ is the only vector that satisfies the associated population first-order conditions, then $b = \beta$ is again the unique solution to the above population minimization problem.

Now form sample counterparts of both $E \left[\phi(X_t,b)\right]$ and $\mathbb{V}(b)$ as functions of $b$. Minimizing a sample counterpart of the above population minimization problem gives rise to a "continuously-updated GMM estimator". See {cite}`HansenHeatonYaron:1996`. The parameter vector and an appropriately scaled minimized objective function have the same limiting distributions as those described in {prf:ref}`rem:old5.2`.[^gmmFootnote]

[^gmmFootnote]: {cite}`ChernozhukovHong:2003` and {cite}`Chen2018` devise and justify simulation-based methods for inference applicable for the continuously-weighted GMM objective function. They do so by adapting insights from simulation-based approaches for Bayesian inferences. Such methods make the statistical analysis of nonlinear moment condition models more tractable.
````
````{prf:example}
Consider the following conditional moment restriction:
```{math}
E \left( {Y_t}'\alpha \mid {\mathfrak A}_{t-\ell} \right) = 0
```
where the random vector $Y_t$ and parameter vector $\alpha$ are both $k \times 1$. We want to know whether there is an $\alpha \ne 0$ that satisfies the conditional moment conditions. Evidently the parameter vector $\alpha$ is only identified up to scale so that the conditional moment restrictions at most identify a one-dimensional family of parameter vectors. In practice, researchers typically achieve identification by normalizing. This can be done, for example, by arbitrarily setting a particular component of $\alpha$ to be unity or else by restricting the norm of $\alpha$ to be unity. If one restricts the norm of $\alpha$ in this way, at best there will be two solutions, say, $\alpha$ and $-\alpha$ with the same norms. Economic interpretations should guide a normalization.

By taking an $r$ dimensional vector $Z_{t-\ell}$ that is in the conditioning information set at date $t-\ell$ and thus ${\mathfrak A}_{t-\ell}$ measurable, we can form an implied unconditional moment condition,
```{math}
E \left( Z_{t-\ell} {Y_t}' \right) \alpha = 0
```
from which we deduce that $\alpha$ must be in the null space of the matrix
```{math}
E \left( Z_{t-\ell} {Y_t}' \right) .
```
To identify a one-dimensional null space of $\alpha$ vectors, it is necessary that the matrix
```{math}
E \left( Z_{t-\ell} {Y_t}' \right)
```
have rank $k$, a restriction that it is straightforward to test.

Two-stage least square imposes a normalization that affects the one-dimensional null space. A one-dimensional null space is also affected when we use a fixed covariance matrix ${\mathbb V}$. In contrast, normalizations imposed in "continuously updated GMM" typically do not affect the null space.
````

<!-- 
## GMM Efficiency Bound

```{math}
\sqrt{N} (b_N - \beta) \approx \textrm{ Normal} \left( 0 , cov(A) \right)
```
where
```{math}
cov(A) = (A' D)^{-1} A'VA (D' A)^{-1}
```
Find a greatest lower bound for $cov(A)$ provided that there exists an $A$ such that $D'A$ is nonsingular.

- Impose $D'A = I$ without loss of generality

- Find ${\widetilde A}$ such that 
```{math}
A' D = A' V {\widetilde A}
```
for all $A$.  

- Construct 
```{math}
A^*  =  {\widetilde A} \left( D' {\widetilde A} \right)^{-1}
```
and note that $D'A^* = I$.

- Show that 
```{math}
cov(A) \ge \left( D' V^{-1} D \right)^{-1} = cov\left( {\widetilde A} \right) = cov \left( A^* \right)
```
## Approximation of Sample Moment Conditions

Use 

```{math}
\frac 1 {\sqrt{N}} \sum_{t=1}^N F(X_t, b_N) 
```

taking account of the estimation of $\beta$

- Recall

```{math}
\frac 1 {\sqrt{N}}  \sum_{t=1}^N  F(X_t , \beta)  
\approx \text{Normal} ( 0 , V) 
```

- 

```{math}
\frac 1 {\sqrt{N}} \sum_{t=1}^N F(X_t, b_N) \approx \left[ I - D (A' D)^{-1} A' \right] \frac 1 {\sqrt{N}}  \sum_{t=1}^N  F(X_t , \beta) 
```

- Observation

```{math}
A' \left[ I - D (A' D)^{-1} A' \right] = 0.
```
## Introduction

The theoretical and applied econometrics literature makes repeated references to generalized method of moments (GMM) estimation as analyzed in {cite}`hansen82`. {cite}`hansen85` lays out a much richer specification of GMM estimators than {cite}`hansen82`, often better matched to the underlying econometric restrictions. {cite}`hansen82` used an extension of an approach of {cite}`sargan58` to represent a finite-dimensional family of GMM estimators and study their properties. {cite}`hansen85` constructed a GMM efficiency bound rich enough to accommodate conditional moment restrictions. Both {cite}`hansen82` and {cite}`hansen85` allow explicitly for temporal dependence.

In this note, we study two optimization problems. The first is the GMM efficiency bound.

This note posits a family of GMM estimators and constructs the greatest lower bound on the statistical efficiency. There is a related efficiency problem that considers a family of maximum likelihood estimators and computes the least informative parameterization that satisfies the underlying econometric restrictions. This is a well known way for computing semiparametric efficiency bounds. See {cite}`newey1990` and {cite}`bkrw1993`.
## Basic set up

Following {cite}`hansen85`, let $Z$ be an index set $\mathcal{Z}$ that is a linear space of $m$-dimensional random vectors used to determine the model implications in a GMM fashion.[^hansen_footnote]

Let $\phi$ be an $m$ dimensional function and suppose the economic model implies:
```{math}
:label: basic_model
E \left[Z \cdot \Phi(Y, \beta) \right] = 0
```
for all $Z \in \mathcal{Z}$ where $\beta$ is an unknown parameter.
Since the parameter vector is $k$-dimensional, for estimation we will use $Z^j$ entries in $\mathcal{Z}$ for $j=1,2,.., k$,
which gives us $k$ equations in $k$ unknowns:
```{math}
E \left[Z^j \cdot \Phi(Y, b) \right] = 0
```
to use in identifying $b = \beta$.

We consider three important special cases.

### Conditional moment restrictions

Suppose that
```{math}
:label: conditional_moment
E\left[ \Phi(Y, \beta) \vert \mathfrak{K} \right] = 0
```
for some conditioning information set $\mathfrak{K}$.  Suppose that $\phi(X, b)$ for $b \in \mathbb{R}^k$ has a finite second moment.  Let $\mathcal{Z}$ contain at least $m$-dimensional random vectors that are bounded and $\mathfrak{K}$ measurable.  Then the family of unconditional moment restrictions of the form {eq}`basic_model` for $Z \in \mathcal{Z}$ imply the conditional moment restrictions.

### Limited number of unconditional moment restrictions

Suppose that $\mathcal{Z}$ consists of linear combinations of $r \ge k$ basis $Z$'s that satisfy {eq}`basic_model`. Stack the $r$ versions of {eq}`basic_model` and write them as:
```{math}
E\left[ F(X, \beta) \right] = 0
```
where $X$ contains $Y$ and the $r$ basis $Z$'s. Now let $\mathcal{Z} = \mathbb{R}^r$. In this case an estimator may be associated with an $r \times k$  selection matrix $A$  or real numbers where each of the $k$ columns are vectors in $\mathcal{Z} = \mathbb{R}^r$. Thus the equations used to estimate are:
```{math}
A' E\left[ F(X, b) \right] = 0
```
which are satisfied for $b = \beta$. 

### Moment matching

Suppose that
```{math}
\Phi(Y, b) = \Psi(Y) - \Gamma(b) 
```
where
```{math}
E\left[\Psi(Y) \right] = \Gamma(\beta).
```
In this formulation, $\Psi(Y)$ defines the moments to be matched and $\Gamma(b)$ gives the model predicted moments as a function of a potential parameter vector $b$. We presume that $r \ge k$. Again $\mathcal{Z} = \mathbb{R}^r$.

[^hansen_footnote]: {cite}`hansen85` allowed the indexes to be random matrices instead of vectors, but this distinction will not be important in what follows.
## Data generation

Next we introduce the data evolution using a measure-preserving transformation $ {\mathbb S}$ where

```{math}
X_t = X \circ {\mathbb S}^t \\
Y_t = Y \circ {\mathbb S}^t \\
Z_t = Z \circ {\mathbb S}^t
```

The resulting processes are stationary. As we will be using time series to approximate expectations, we should think of this expectation as conditioned on invariant events along with the parameter $\beta$.
Mean-value approximation

It suffices for us to consider estimators that satisfy 
```{math}
A' {\frac 1 {\sqrt{N}}} \sum_{t=1}^N F(X_t, b_N) = 0.
```
Let $D$ be the $r$ by $k$ vector of derivatives:
```{math}
D = E\left[ \frac {\partial F}{\partial b'}(x, \beta) \right].
```
Then an application of the Mean-Value Theorem along with the Law of Large Numbers for stationary processes
informs us that 
```{math}
0 \approx A' {\frac 1 {\sqrt{N}}} \sum_{t=1}^N F(X_t, \beta) + A' D \sqrt{N} (b_N - \beta)
```
or
```{math}
\sqrt{N}(b_N - \beta)  \approx - (A' D)^{-1} A' {\frac 1 {\sqrt{N}}} \sum_{t=1}^N F(X_t, \beta).
```
provided that $A'D$ is nonsingular. This approximation was established formally in {cite}`hansen82`.


````{prf:theorem}
:label: assume:identify
There is a matrix $A$ such that $A'D$ is nonsingular.
````

This presumes that the model is "locally identified".
## Martingale approximation

Following {cite}`gordin` and {cite}`hansen85` map the stochastic process $\{F(X_t, \beta)\}$ into a sequence of martingale differences $\{H_t\}$ with an equivalent central limit approximation. Let

```{math}
U_t = F(X_t, \beta)
```

Construct

```{math}
G_t = \sum_{j=0}^\infty E \left[U_{t+j}\vert {\mathfrak F}_{t} \right]
```

where we presume the right-hand side converges in mean-square. Form the innovation:

```{math}
H_{t+1}  = G_{t+1} - E\left(G_{t+1} \vert {\mathfrak F}_{t} \right).  
```

Note that

```{math}
U_t =  G_t - E\left(G_{t+1} \vert {\mathfrak F}_t \right) = H_{t+1} + G_t - G_{t+1}.
```

Then

```{math}
{\frac 1 {\sqrt{N}}} \sum_{t=1}^N U_t \approx {\frac 1 {\sqrt{N}}} \sum_{t=1}^N H_{t+1} \Rightarrow \textrm{Normal}\left(0, 
V \right)
```

where

```{math}
V = E\left( H_{t+1} {H_{t+1}}' \right).  
```


````{prf:theorem}
The covariance matrix $V$ is nonsingular.  
````
## Finite-dimensional GMM efficiency bound

Think of the selection matrix $A$ as a way to index alternative GMM estimators. Given the Central Limit approximation, construct an asymptotic covariance matrix for $\{A \frac{1}{\sqrt{N}} \sum_{t=1}^N U_t : N=1,2,... \}$ as
```{math}
< A | A > = A' V A
```
Now represent
```{math}
A' D = A' V A^* = < A | A^* >
```
for all $A$ where $A^* = V^{-1} D$. The Mean Value Approximation implies that the asymptotic covariance matrix for the associated GMM estimator is:
```{math}
cov(A) = (< A | A^* >)^{-1} < A | A > (< A^* | A >)^{-1}
```

````{prf:theorem}
The GMM efficiency bound is:

```{math}
cov(A) \ge (A^* | A^*)^{-1} = cov(A^*) = (D' V^{-1} D)^{-1} \overset{\text{def}}{=} inf
```

A GMM counterpart to information is given by the respective entries of $(A^* | A^*)$.

**Proof.** We establish the proposition in two steps. First let $B$ be a nonsingular $k \times k$ matrix. Note that 

```{math}
(A B | A^*) = B'(A' | A^*)
```

and 

```{math}
\left(B'(A' | A^*)\right)^{-1} = \left((A' | A^*)\right)^{-1} \left(B' \right)^{-1}
```

It follows that $cov(AB) = cov(A)$. This should not be surprising as premultiplying an equation system by a nonsingular matrix $B'$ does not alter the solution to the equation system. Without loss of efficiency, we may impose that
$(A | A^*) = I$.  

Next, let $\widetilde{A} = A^* (A^* | A^*)^{-1}$. It may be shown that if $(A | A^*) = I$,

```{math}
(A - \widetilde{A} | A - \widetilde{A}) = (A | A) -  (\widetilde{A} | \widetilde{A}) = (A | A) - inf.
```

This verifies the efficiency bound the left-hand matrix is positive semidefinite.
````
## Extending this bound to a larger set of moment conditions

Following {cite}`hansen85`, we consider a more general index set $\mathcal{Z}$. We use the same two approximations as for the finite-dimensional case.

### Martingale approximation

For a $Z \in \mathcal{Z}$, form the scalar $R_t$

```{math}
R = Z \cdot \phi(Y, \beta),
```

and the infinite sum

```{math}
M = \sum_{j=0}^\infty E\left( R_{j+1} \vert \mathfrak{F}_1\right)  - \sum_{j=0}^\infty E\left( R_{j+1} \vert \mathfrak{F}_0\right)
```

used in the martingale approximation. By construction:

```{math}
E\left(M \vert \mathfrak{F}_0 \right) = 0
```

and thus $\{ \sum_{t=1}^{N}  M_{t} \}$ is a martingale.

Now define linear mapping:

```{math}
\mathbb{M}(Z) = M
```

taking the random vector $Z$ into the scalar random variable $M$, and let $\mathcal{N}$ be the linear space of random variables:

```{math}
\mathcal{N}^o = \{ \mathbb{M}(Z) : Z \in \mathcal{Z} \}.
```

Construct $\mathcal{N}$ as the mean square closure of $\mathcal{N}$.

### Including the mean value approximation

Define:

```{math}
D^j = E\left[ \frac {\partial \phi}{\partial b_j}(x, \beta) \right].
```

where $b_j$ is the $j^{th}$ coordinate of the $k$-dimensional vector $B$. Introduce the linear functionals on $\mathcal{Z}$

```{math}
\mathbb{L}^j(Z) = E(Z \cdot D^j)
```


````{prf:theorem}
:label: assumption-mean-value-approximation

For any sequence $\{Z^i : i=1,2,... \}$ for which $\{ \mathbb{M}(Z^i) : i=1,2, ... \}$ converges in mean square to zero, $\{ E \left( Z^i \cdot D^j\right) : i=1,2,...\}$ converges to zero for each $j=1,2,...,k$.
````

Under this assumption, by the Riesz Representation Theorem, there exist $N^j$ in $\mathcal{N}$ such that

```{math}
E Z \cdot D^j = E\left[\mathbb{M}(Z) N^j\right]
```

for all $Z$ in $\mathcal{Z}$ and $j=1,2,...,k$.

A GMM estimator is constructed by choosing $Z^i$ for $i = 1,2,...,k$. Adapting and modifying notation a bit, the covariance matrix for this estimator is:

```{math}
\begin{align*}
cov(Z^1,Z^2,.., Z^k)  = &\left[E\left( \begin{bmatrix} {\mathbb M}(Z^1) \cr {\mathbb M}(Z^2) \cr ... \cr {\mathbb M}(Z^k) \end{bmatrix}
\begin{bmatrix} N^1 & N^2 & .... &N^k \end{bmatrix} \right)\right]^{-1}  \cr & E\left( \begin{bmatrix} {\mathbb M}(Z^1) \cr 
{\mathbb M}(Z^2) \cr ... \cr {\mathbb M}(Z^k) \end{bmatrix} 
\begin{bmatrix} {\mathbb M}(Z^1) & 
{\mathbb M}(Z^2) & ... & {\mathbb M}(Z^k)  \end{bmatrix} \right) \cr & \left[E\left( \begin{bmatrix} N^1 \cr N^2 \cr ....  \cr N^k \end{bmatrix} 
\begin{bmatrix} {\mathbb M}(Z^1) & {\mathbb M}(Z^2) & ... & {\mathbb M}(Z^k) \end{bmatrix}
 \right)\right]^{-1}
\end{align*}
```

````{prf:proposition}
```{math}
cov(Z^1,Z^2,.., Z^k) \ge \left[ E\left( \begin{bmatrix} N^1 \\ N^2 \\ .... \\ N^k \end{bmatrix} \begin{bmatrix} N^1 & N^2 & .... & N^k \end{bmatrix} \right)\right]^{-1} \overset{\text{def}}{=} inf
```
````

:::{proof}
The bound follows from essentially the same argument as in the finite dimensional case. It may be shown that the bound is sharp by using mean square limiting arguments.
:::

---

With this construction, we form

```{math}
{\mathbb M}(z) = h
```

where $E(h_{t+1} | {\mathcal F}_t ) = 0$.

With this construction:

```{math}
\frac 1 {\sqrt{N}} \sum_{t=1}^N z_t \cdot f(y_t, \beta) \approx \frac 1 {\sqrt{N}} \sum_{t=1}^N {\mathbb M}(z)_t ,
```

where

```{math}
\frac 1 {\sqrt{N}} \sum_{t=1}^N {\mathbb M}(z)_t \Rightarrow N \left( 0, {\mathbb M}(z)^2 \right),
```

and $\Rightarrow$ denotes weak convergence.

Define the inner product.

```{math}
< z | z^* > = E\left[ {\mathbb M}(z) {\mathbb M}(z^*) \right],
```

Let

```{math}
{\mathcal U}^o = \left\{ u = {\mathbb M}(z) : z \in {\mathcal Z} \right\}
```

and ${\mathcal U}$ the mean square closure of ${\mathcal U}^o$.

Finally define

```{math}
{\mathbb M}(Z) = \begin{bmatrix} {\mathbb M}(z^1) \\ {\mathbb M}(z^2) \\ ... \\ {\mathbb M}(z^k) \end{bmatrix}
```
## Combining the approximations

For a random matrix $Z$ with columns in ${\mathcal Z}$, the corresponding GMM estimator satisfies:
```{math}
\sqrt{N}(b_N - \beta) \approx  - \left[E\left(Z' D \right)\right]^{-1} {\frac 1 {\sqrt{N}}} \sum_{t=1}^N {\mathbb M}(Z)_t
```
We define the corresponding asymptotic covariance matrix as:
```{math}
 {\mathbb A}(Z) = \left[E\left(Z' D \right)\right]^{-1}E\left[{\mathbb M}(Z){\mathbb M}(Z)'\right] \left[E\left(D' Z \right)\right]^{-1}
```
where we have denoted explicitly the dependence of the asymptotic covariance matrix on the choice of $Z$. The function
${\mathbb A}$ maps $Z$'s into positive semi-definite matrices.

Prior to computing the GMM efficiency bound, we obtain an alternative representation of $EZ'D$.
Consider first the $k$ linear functionals of ${\mathcal Z}$:
```{math}
E z \cdot d^j
```
where $d^j$ is the $j^{th}$ column of $D$:
```{math}
d^j = {\frac {\partial f}{\partial b^j}}
```
where $b^j$ is the $j^{th}$ coordinate of $b$.
We presume that If ${\mathbb M}(z) =0$, then $Ez \cdot d^j = 0$ for $j=1,2,...,k$.  We strengthen this by assuming:


````{prf:assumption}
For any sequence $\{z^i : i=1,2,... \}$ for which $\{ {\mathbb M}(z^i) : i=1,2, ... \}$ converges in mean square to zero, $\{ E z^i \cdot d^j : i=1,2,...\}$ converges to zero for each $j=1,2,...,k$.
````

Under this assumption, by the Riesz Representation Theorem, there exist $u^j$ in ${\mathcal U}$ such that
```{math}
E z \cdot d^j = E\left[{\mathbb M}(z) {\tilde u}^j\right]
```
for all $z$ in $\mathcal{Z}$ and $j=1,2,...,k$.  Construct:
```{math}
{\widetilde U} = \begin{bmatrix} {\tilde u}^1 \\ {\tilde u}^2 \\ ... \\ {\tilde u}^k \end{bmatrix}.
```
Then
```{math}
E\left(Z'D\right) = E\left[ {\mathbb M}(Z) {\widetilde U}' \right],
```
In light of Assumption {prf:ref}`assume:identify`, $E \left( {\widetilde U} {\widetilde U}' \right)$ is nonsingular.  
Moreover, 
```{math}
 {\mathbb A}(Z) =  \left(E\left[U {\widetilde U}'\right]\right)^{-1} E \left( U U' \right) \left[E\left({\widetilde U} U' \right)\right]^{-1}
```
for $U = {\mathbb M}(Z)`.
## GMM efficiency bound

Prior to establishing the efficiency bound, we study the function

```{math}
{\mathbb B}(U) = \left[E\left(U {\widetilde U}'\right)\right]^{-1} E \left( U U' \right) \left[E\left({\widetilde U} U' \right)\right]^{-1}
```
where the entries of $U$ are in ${\mathcal U}$. Notice that ${\mathbb B}$ maps random vectors $U$ with entries in 
${\mathcal U}$ into positive semidefinite matrices.  

Let $L$ be a $k$-dimensional nonsingular matrix of real numbers. Then

```{math}
{\mathbb B}(LU) = {\mathbb B}(U)   
```
In light of this, when we seek to find a greatest lower bound for ${\mathbb B}$, it suffices to consider only $U$'s such that 

```{math}
E\left(U {\widetilde U}'\right) = I
```
Form

```{math}
{\widehat U} = \left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}{\widetilde U}.
```
Then

```{math}
{\mathbb B}\left({\widehat U} \right) = E\left( {\widehat U}{\widehat U}' \right) = \left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}.
```

For $U$ such that $E\left(U {\widetilde U}'\right) = I$, compute

```{math}
E\left[ \left(U - {\widehat U}\right)\left(U - {\widehat U}\right)' \right] = {\mathbb B}(U) - 2\left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}
+ {\mathbb B}({\widehat U}) = {\mathbb B}(U) -  {\mathbb B}({\widehat U})
```
where the left-hand side matrix is positive definite. Thus

```{math}
 {\mathbb B}(U) \ge  \left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}
 ```
where the inequality $\ge$ implies that the matrix on the left-hand side minus the matrix on the right-hand side is positive semidefinite. It follows that

```{math}
{\mathbb A}(Z) \ge \left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}
```
for $Z$ with columns in ${\mathcal Z}$. This bound can be attained since we may construct a sequence $\left\{ Z^i : i=1,2,...\right\}$
such that $\left\{ {\mathbb M}\left(Z^i\right) : i=1,2,...\right\}$ converges in mean square to ${\widetilde U}$. Thus 
$\left[E\left({\widetilde U} {\widetilde U}'\right)\right]^{-1}$
 is a sharp efficiency bound for the family of GMM estimators.
## Efficiency result from {cite}`hansen82`

Let 

```{math}
E\left[ f(x, \beta) \right] = 0
```
be a vector of unconditional moment conditions with $r \ge k$. In this case $z$ is an $r$-dimensional vector of real numbers and can be brought out of the expectation. This unconditional moment condition could be derived from a conditional moment condition

```{math}
E\left[\phi(y, \beta) \vert {\mathfrak G} \right] = 0
```
where $\phi$ has $q$ coordinates. Construct an $r \times q$ matrix $X$ with entries that are ${\mathfrak G}$ measurable implying that:

```{math}
E\left[ X \phi(y, \beta) \right] = 0. 
```
In this case we let 

```{math}
:label: junk1
f(x,b) = X \phi(y, b)
```
where the vector $x$ contains the entries of $y$ and $X$.

Form

```{math}
H_t = \sum_{j=0}^\infty \left( E\left[ f(x_{t+j}, \beta) \vert {\mathcal F}_t \right] - E\left[f(x_{t+j}, \beta) \vert {\mathcal F}_{t-1} \right]\right)
```
used for a martingale approximation for 

```{math}
\frac 1 {\sqrt{N}} \sum_{t=1}^N f(x_t, \beta).
```
Let 

```{math}
V = E\left( H H' \right),
```
which we take to be nonsingular. Then

```{math}
{\mathcal M}(Z) = Z'H.
```
Moreover,  

```{math}
{\widetilde U} = {\widetilde Z}' H
```
satisfies:

```{math}
Z' {\overline D} = Z' V {\widetilde Z}.
```
where ${\overline D} = E\left(D\right)$ for all $r \times k$ matrices $Z$. Since this applies for all $r \times k$ matrices $Z$,

```{math}
Z^* = V^{-1} {\overline D} 
```
and the efficiency bound is $\left({\overline D}' V^{-1}{\overline D} \right)^{-1}$.
## Conditional moment restriction

In construction {eq}`junk1`, we had great flexibility in the construction of the matrix $X$. We now want to exploit that flexibility to improve the asymptotic efficiency of the resulting estimator. With this in mind, we now set 

```{math}
f(x, b) = \phi(y, b),
```

and hence we consider a conditional moment restriction of the following form:

```{math}
E\left[   f(y, \beta) \vert {\mathfrak G} \right] = 0.  
```

This is equivalent to the supposition that

```{math}
E\left[  z \cdot  f(y, \beta)  \right] =0 
```

where the entries of the $z$'s are functions that are measurable with respect to ${\mathfrak G}$.  

We suppose either that the data generation is iid or more generally that the conditioning information in ${\mathfrak G}$ includes the past time series so that

```{math}
E\left[z_t \cdot  f(y_t, \beta) \vert {\mathfrak F}_{t-1} \right] = 0.
```

This simplification implies that 

```{math}
{\mathbb M}(Z) = Z f(x, \beta) 
```

We proceed by solving:

```{math}
z ' E\left(d^j \vert {\mathcal G}  \right) = z' E\left[ f(y, \beta) f(y, \beta)' \vert {\mathcal G} \right]  {\tilde z}^{j},
```

for all admissible $z$. Solving for ${\tilde z}^j$ gives: 

```{math}
{\tilde z}^{j} = \left( E\left[ f(y, \beta) f(y, \beta)' \vert {\mathcal G} \right] \right)^{-1} E\left(d^j \vert {\mathcal G}\right).
```

Thus 

```{math}
{\tilde u}^j = {\tilde z}^j \cdot f(y, \beta) .
```

The resulting $GMM$ efficiency bound is:

```{math}
\left(E \left[  E\left(D' \vert {\mathcal G}\right)  \left( E\left[ f(y, \beta) f(y, \beta)' \vert {\mathcal G} \right] \right)^{-1} E\left(D \vert {\mathcal G}\right) \right] \right)^{-1}
```
## Relations to maximum likelihood

Suppose the data generation is iid. Consider first the case where $z$ is restricted to be a $k$-dimensional vector of numbers and thus the moment conditions of interest are:
```{math}
E\left[ f(x, \beta) \right] = 0
```
Solve the problem:

````{prf:problem}
:label: expotilt
\[
\min_{m, m\ge0, Em=1} E(m \log m ) 
\]
subject to
\[
E[ m( f(x,b) ] = 0.
\]
````
Here $m$ is a relative density for $x$ and it depends on $x$ and a hypothetical parameter vector $b$. The objective is to minimize what is called  *relative entropy* subject to the moment conditions. 
This is a population version of the relative entropy problem proposed by {cite}`kitamurastutzer97`. It is interesting but intuitive that the solution to this problem provides the least informative parameterization of the parameter of interest. We solve this problem for an arbitrary $b$ perhaps restricted to be in a neighborhood of $\beta$.

The solution to Problem {prf:ref}`expotilt` is well known to entail exponential tilting
```{math}
m^*(x, b ) = {\frac {\exp[ \lambda(\theta) \cdot f(x, b) ] }{E \left( \exp[ \lambda(\theta) \cdot f(x,b) ] \right)}} .
```
where $\lambda(b)$ is chosen to satisfy:
```{math}
\max_{\lambda} - \log E \left( \exp \left[ \lambda \cdot f(x, b) \right] \right)
```
Maximizing over $\lambda$ ensures that the moment conditions are satisfied and $\lambda(\beta) = 0$.  

We take $m^*(x, b)$ to be density relative to the true data generating process. Suppose we construct a log-likelihood with this density and use to estimate $\beta$.     
The score with respect to $b$ is
```{math}
S(x)  = \left[{\frac {\partial \lambda}{\partial b}}(\beta) \right]'  f(x,\beta). 
```
To compute ${\frac {\partial  \lambda}{\partial b}}$, we differentiate inside the expectation
```{math}
E\left[ m^*(x, b) f(x,b) \right] = 0
```
to show that 
```{math}
E\left[ S(x) f(x, \beta)'\right] + {\overline D} '  = 0.  
```
Substituting for $S$, we see that 
```{math}
 \left[{\frac {\partial \lambda }{\partial b}}(\beta)\right]'V = - {\overline D}', 
```
or
```{math}
\left[{\frac {\partial \lambda }{\partial b}}(\beta)\right] = V^{-1} {\overline D}.
```
Thus
```{math}
S(x) = - {\overline D}' V^{-1} f(x, \beta) = - {\widetilde U}(x).
```
The efficiency of this likelihood estimator is:
```{math}
\left[E \left(S S' \right)\right]^{-1} = \left[E \left( {\widetilde U}{\widetilde U}' \right)  \right]^{-1}
```
which is the same as the GMM efficiency bound.  

What do we make of this calculation? Consider parameterizations of the density for $x$  for which the moment conditions are satisfied. The efficiency of such estimators should exceed that of all of the GMM estimators. Since we found a parameterization for which ML efficiency is the same as the GMM efficiency bound we know the GMM lower bound coincides with the ML upper bound.  

We perform the analogous calculation in terms of conditional moment restriction to parameterize the conditional density for $x$ given ${\mathfrak G}$. We construct $\lambda(b)$ as a function of the conditioning information set and we 
solve:
```{math}
 \left[{\frac {\partial \lambda }{\partial b}}(\beta)\right]'E\left[ f(x, \beta) f(x,\beta)'  \vert {\mathfrak G} \right]  = - E\left(D' \vert {\mathfrak G} \right) 
```
Thus
```{math}
S(x) = - E\left(D' \vert {\mathfrak G} \right) \left(E\left[ f(x, \beta) f(x,\beta)'  \vert {\mathfrak G} \right]\right)^{-1}
f(x, \beta) = - {\widetilde U}(x).  
```
Which again shows the connection between ML and GMM efficiency.  
Such a calculation provides an informal derivation of {cite}`chamberlain87`'s semiparametric bound for conditional moment restrictions.  



We now consider a family of estimators that applies this approach for alternative choices of $z^1, z^2, ..., z^k$.  
From {cite}`hansen85`, the GMM efficiency bound is
```{math}
\left( E \left[ {\rm vec}(g^i){\rm vec}(g^i)' \right]  \right)^{-1} 
```
By a partitioned inverse formula, the reciprocal of the GMM efficiency bound for the first entry of $\beta$ is given by the inverse of the regression error variance obtained by a population least squares regression of ${g}^1$ onto $g^2$, ..., ${g}^k$. Consequently we compute
```{math}
:label: solution
u = {g}^1 - Proj \left({g}^1 \vert {g}^2, {g}^3, ..., {g}^k \right)
```
where $Proj(\cdot)$ is the least squares projection operator. Then ${\frac 1 {E\left( {u}^2 \right)}}$ is the GMM efficiency bound for the first entry of the parameter vector.   


````{prf:problem}
\[
\inf_{f^1, f^2, f^k \in {\mathbb F}^o}   \begin{bmatrix} 1 & 0 & ... & 0 \end{bmatrix} 
 E \left[{\rm vec}(f^i) {\rm vec}(f^i)' \right]
\begin{bmatrix} 1 \cr 0 \cr ... \cr 0 \end{bmatrix}
\]
subject to $E\left[{\rm vec}(f^i) {\rm vec}(g^i)'\right] = I$.  
````
The identity matrix constraint is imposed for convenience to simplify the objective. Imposing this constraint does not alter the implied efficiency bound because we may always take linear combinations of the $z^i$ to satisfy the constraint without altering the implied GMM estimator.  
The infimum is given by ${\frac 1 {E u^2}}$. By replacing ${\mathbb F}^o$ with its closure ${\mathbb F}$, the infimum is attained by
```{math}
{\rm vec}(f^i) =  \left(E \left[{\rm vec}(g^i) {\rm vec}(g^i)' \right]\right)^{-1}  {\rm vec}(g^i)
```
provided that ${\rm vec}(g^i)$ has a nonsingular second moment matrix. This solution remains the same no matter which linear combination of the coefficients we feature in the minimization, although the implied infimum will of course be different.  

To summarize, i) we represent the partial derivatives of the moment condition for alternative choices of $z$ as a bounded linear functional on ${\mathbb F}$, ii) we represent these functionals using the Riesz Representation Theorem, and iii) we run a population least squares regression to obtain the efficiency bound for each parameter.
## Efficiency result from {cite}`hansen82`

The $z$'s are vectors of real numbers can be brought out of the expectation. Then

```{math}
z' E d^i = z' V \tilde{z}^i
```

for all $z$ where

```{math}
V = E \left[f(y, \beta) f(y, \beta)'\right].
```

Thus, provided that $V$ is nonsingular:

```{math}
g^i = \tilde{z}^i \cdot f(y, \beta)
```

where

```{math}
\tilde{z}^i = V^{-1} E d^i.
```

This gives the GMM efficiency bound computed in {cite}`hansen82`:

```{math}
\left[ E (d') V^{-1} E d \right]^{-1}.
```
## Conditional moment restriction

Suppose the entries of the $z$'s are functions that are measurable with respect to $\mathcal{B}$ where 
```{math}
E\left[ f(y, \beta) \vert \mathcal{B} \right].
```
We proceed as in the previous case by solving:
```{math}
z ' E\left(d^i \vert \mathcal{B}  \right) = z' E\left[ f(y, \beta) f(y, \beta)' \vert \mathcal{B} \right] \tilde{z}^i,
``` 
for all $z$ or 
```{math}
\tilde{z}^i = \left( E\left[ f(y, \beta) f(y, \beta)' \vert \mathcal{B} \right] \right)^{-1} E\left(d^i \vert \mathcal{B}\right).
```
Thus 
```{math}
g^i = \tilde{z}^i \cdot f(y, \beta) .
```
## Dual efficiency problem

We next compute the least informative parameterization of the model that satisfies moment conditions {eq}`basicmodel`.  
Consider a density over the data parameterized by a scalar $\sf{r}$ where $\sf{r} = \theta $ is the value that gives the true density.    Let $s$ be the corresponding score random variable and
let $b(\sf{r})$ be the implied choice of the parameter vector $\beta$.  In particular,
$b(\theta) = \beta$.  Define:
```{math}
\alpha \cdot \frac {\partial b}{\partial \sf{r}} (\theta)
```
Then by standard likelihood calculations, the Fisher information for the estimator of $\theta$ is
$E s^2$ and the implied asymptotic covariance matrix for the estimator of $\beta$ is
```{math}
{\frac 1 {E s^2}} \alpha \alpha'
```

Moment restrictions implicitly restrict the score.  
An integration-by-parts argument implies 
```{math}
E \left[z \cdot f(y,\beta) s \right] =  E f s  = - \sum_{i=1}^k \alpha_i D_i(f)  = - \sum_{i=1}^k \alpha_i  E f g^i
```
for any $f \in {\mathbb F}$.    Thus {eq}`basicmodel` restricts the score $s$ to satisfy:
```{math}
:label: scorerestrict
Proj( s \vert {\mathbb F} ) =  - \sum_{i=1}^k \alpha_i  g^i .
```
Another way to state the restriction on the score is to let 
${\mathbb G}$ be the linear space generated by $g^i$ for $i=1,2,...k$ and to require that
```{math}
:label: restrict
Proj( s \vert {\mathbb F} ) =  Proj \left( s \vert {\mathbb G} \right).
```
In what follows,
suppose now that dual representation of the model as a restriction {eq}`restrict` on the scores.  We have shown that scores should satisfy this restriction, but for the moment we have not explored the converse. We have not taken any such candidate score and produced an implied parameterization
of the model subject to the moment conditions.  We will have more to say about this later.  

{cite}`chamberlain87` established an efficiency bound for conditional moment restrictions
using multinomial approximation and GMM estimation.  Here we 
study of semiparametric efficiency by computing the least information parameterizations.  
We consider this approach in what follows.   Given the fully parameterized model
of the probability distribution, the reciprocal of the  implied variance for the first entry of $\beta$ is 
```{math}
{\frac {(\alpha_1)^2}{ Es^2}} 
```
In light of this formula,  the least informative scalar parameterization has a score that solves:


````{prf:problem}
:label: 
\[
\min_{s \in {\mathbb L}^2} E s^2
\]
subject to $Es = 0$ and:
\[
Proj (s \vert {\mathbb F} ) =  {g}^1 - \sum_{i=2}^k   \alpha_i g^i.
\]
````

Since $Proj (s \vert {\mathbb F} )$ has a smaller variance than $s$, it suffices to
consider only $s$ in $ {\mathbb F}$, which by construction have projection errors equal to zero.  
Given the restriction, the solution for $s$ is
```{math}
v = { g}^1 - Proj \left({g}^1 \vert {g}^2, {g}^3, ..., {g}^k \right).
```
where the $\alpha_i$'s are the population regression coefficients.  
Notice that ${u} = {v}$ used to represent the GMM efficiency bound.  

Let me now take inventory.  We have derived a necessary condition for the restrictions on 
scalar score random variables.  We not yet taken any such candidate score and produced an implied parameterization
of the model subject to the moment conditions.  We use the necessary condition to restrict the scores, and in so doing provide 
a direct connection to the GMM efficiency bound and a semiparametric counterpart.  With more specificity on the problem, we may sidestep the need to map potential scores into parameterizations that satisfy the moment conditions.
## Parameterizing the Model

We consider two cases. Rather than approximating the entire family of scores, we show how to construct least informative ones by solving minimum relative entropy problems.

### Unconditional moment conditions

Consider first the case where $z$ is a $k$-dimensional vector of numbers. Depict the parameter of interest as a function $\beta + \theta$. Solve the problem:


````{prf:problem}
:label: expotilt

```{math}
\min_{m, m\ge0, Em=1} E(m \log m ) \text{ subject to } E[ m( f(y, \beta + \theta) ) ] = 0. 
```
````

Here $m$ is a relative density for $y$ and it depends on $y$ and a hypothetical parameter vector $b$. The objective is to minimize what is called *relative entropy* subject to the moment conditions. This is a population version of the relative entropy problem proposed by {cite}`kitamurastutzer97`. It is interesting but intuitive that the solution to this problem provides the least informative parameterization of the parameter of interest.

The solution to Problem {prf:ref}`expotilt` is well known to entail exponential tilting

```{math}
m^*(y, \theta ) = \frac {\exp[ \lambda(\theta) \cdot f(y + \theta) ] }{E \left( \exp[ \lambda(\theta) \cdot f(y, \beta + \theta)] \right)} .
```

where $\lambda(\theta)$ is chosen to satisfy:

```{math}
\max_{\lambda} - \log E \left( \exp \left[ \lambda \cdot f(y, \beta + \theta) \right] \right)
```

Maximizing over $\lambda$ ensures that the moment conditions are satisfied and $\lambda(0) = 0$. The score with respect to $\theta$ is

```{math}
s(y) = \left({\frac {\partial \lambda}{\partial \theta}}\vert_{\theta = 0} \right) f(y,\beta). 
```

To compute $\frac {d \lambda^i}{d \theta}$, we differentiate inside the expectation to show that

```{math}
E\left[ s(y) f(y, \beta)\right] + E d^i = 0.  
```

Substituting for $s^i$, we find

```{math}
 \left[\frac {d \lambda^i}{d \theta}(0)\right]'V = - E(d^i), 
```

or

```{math}
\frac {d \lambda^i}{d \theta}(0) = V^{-1} E(d^i).
```

Thus

```{math}
s^i(y) = - (E d^i)' V^{-1} f(y, \beta) = - g^i (y).  
```

By reproducing this argument for alternative $i$ and look across linear combinations of parameterizations and scores. The least-informative score problem is indeed among the linear combination of such scores. We minimize the variance (Fisher information) to find the least informative such parameterization for say the first entry of $\beta$. This is given by running a regression of $g^1$ onto $g^2$, ... $g^k$ and computing regressions error variance. Thus the solution to the least-informative score problem is indeed among the scores for this problem.

Consider next the case of a conditional moment restriction. We may solve the same problem but a conditional version of it. In the score formula we form conditional expectations and reproduce the comparable result. In particular,

```{math}
s^i(y \vert \mathcal{B}) = - \left[E\left(d^i \vert \mathcal{B} \right)\right]'   \left( E\left[ f(y, \beta) f(y, \beta)' \vert \mathcal{B} \right] \right)^{-1}f(y, \beta) = - g^i 
```

{cite}`backbrown92` previously proposed a different parameterization that gives this same score. {cite}`qinlawless94` (Theorem 3) establish semiparametric efficiency for the empirical likelihood estimator, but under their regularity conditions the first-order efficiency of GMM and empirical likelihood coincide.[^compactSupport]

[^compactSupport]: Using a population counterpart to an empirical likelihood approach to construct a parameterization of the moment conditions as we have one here is problematic without additional compact support conditions on the data distribution.
## Time Series Extension

To produce a time series counterpart requires more effort. Score processes for parametric likelihoods are martingales. In a time series setting {cite}`hansen82` and {cite}`hansen85` both use martingale approximations by applying the central limit theory of {cite}`billingsley` and {cite}`gordin`. Specifically, {cite}`hansen85` maps the stochastic process $\{ z_t \cdot \phi(y_t, \beta)\}$ into sequence of martingale differences $\{f_t \}$ with an equivalent central limit approximation. The inner product and the $\mathbb{L}^2$ are constructed using the corresponding martingale difference sequence. With this construction to build on, we may extend this dual relation to richer time series setting since score processes increments are themselves martingale differences. Thus the martingale approximations used to represent GMM efficiency can be linked to the score processes for time series parametric ML estimators.

Formally, let $\mathbb{S}$ denote a one-to-one, measure-preserving and ergodic transformation with a measurable inverse. Construct processes:
```{math}
z_t = z \circ \mathbb{S}^t \\
y_t = y \circ \mathbb{S}^t. 
```
Construct a filtration $\{ \mathcal{F}_t \}$ where $z_{t-j}$ for all possible $z$ and $j \ge 0$ and $y_{t-j}$ for $j \ge 0$. Construct
```{math}
f_t  = \sum_{j=0}^\infty \mathbb{E} \left[ z_{t+j} \cdot \phi(y_{t+j}, \beta)\vert \mathcal{F}_t \right] - \sum_{j=0}^\infty \mathbb{E} \left[ z_{t+j} \cdot \phi(y_{t+j}, \beta)\vert \mathcal{F}_{t-1} \right]
```
where we presume that 
```{math}
\sum_{j=0}^\infty \mathbb{E} \left[ z_{t+j} \cdot \phi(y_{t+j}, \beta)\vert \mathcal{F}_t \right]
```
converges in mean-square. With this construction, we form
```{math}
\mathbb{M}(z) = f
```
where $\mathbb{E}(f_{t+1} \vert \mathcal{F}_t ) = 0$. We define
```{math}
< z \vert z^* > = \mathbb{E}(f f^*)
```
and apply the Riesz Representation Theorem as before. In the special case in which 
```{math}
\mathbb{M}(z) = z \cdot \phi(y, \beta) 
```
we may proceed as before with conditional moment restrictions by directly constructing least informative scores. But in more general circumstances the Dual efficiency problem comes into play justifying the construction of infinite-dimensional families of parameterizations along with their mean square limits.
## Infinite Dimensional Parameter Space

The econometrics literature has studied semiparametric efficiency in more general circumstances. See {cite}`chamberlain92` and {cite}`aichen03`. Here I merely sketch how the extension works.  
Suppose that 

```{math}
\phi(y, \beta) = \phi(y, \beta_1, \beta_2)
```

where $\beta_2$ is now infinite dimensional. It is now most convenient to work with the counterpart to the dual optimization problem. We need to define a differentiable mapping from $\mathsf{R}$ to $b_1$ and $b_2$ where $b_2$ resides in some appropriate function space.  We represent

```{math}
Efs = - \alpha_1 Ef_{g}^1 -  Ef_{g}^\infty
```

where $g^\infty$ is the counterpart to $\sum_{i=2}^k \alpha_i g^i$ and contained in $\mathbb{F}$. We compute this entity by applying a functional counterpart to the chain rule. Changing the parameterization of the model will alter $\alpha_1$ and $g^\infty$. We form a corresponding space $\mathbb{G}$ and ask that the score $s$ satisfy 

```{math}
Proj(s \vert \mathbb{F}) = Proj \left( s \vert \mathbb{G} \right).  
```

---

GMMbib

---
 -->
