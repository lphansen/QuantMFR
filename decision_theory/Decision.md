# Risk, Ambiguity, and Misspecification: Decision Theory, Robust Control, and Statistics


*Authors: Lars Peter Hansen[^lhansen] and Thomas J. Sargent[^tsargent]*

[^lhansen]: University of Chicago. Email: lhansen@uchicago.edu
[^tsargent]: New York University. Email: thomas.sargent@nyu.edu

---

## Abstract

What are "deep uncertainties" and how should their presence influence prudent decisions? To address these questions, we bring ideas from robust control theory into statistical decision theory. Decision theory has its origins in axiomatic formulations by von Neumann and Morgenstern, Wald, and Savage. After Savage, decision theorists constructed axioms that formalize a notion of ambiguity aversion. Meanwhile, control theorists constructed decision rules that are robust to some model misspecifications. We reinterpret axiomatic foundations of decision theories to express ambiguity about a prior over a family of models along with concerns about misspecifications of the corresponding likelihood functions.

*Keywords: deep uncertainty, ambiguity, misspecification, variational preferences, statistical divergence, relative entropy*

*JEL Codes: C10, C14, C18*

---
(sec:intro)=
## Introduction 

Climate scientists confront "deep uncertainties."[^deepuncertainties] Practicing econometricians also often struggle with uncertainty about their statistical models, but usually with scant guidance from significant advances in decision theory made after {cite}`wald1947essentially,wald1949statistical,Wald_statist_dec_book`, {cite}`Savage:1954`, and {cite}`ellsberg1961risk` because so much recent formal theory of decision making under uncertainty in economics is not cast explicitly in terms of the likelihoods and priors that are the foundations of statistics and econometrics.[^econometricians] Likelihoods are probability distributions conditioned on parameters while priors describe a decision maker's subjective belief about parameters.[^likelihoodspriors] By distinguishing roles played by likelihood functions and subjective priors over parameters, this paper aims to bring contributions to decision theory after Wald and Savage into closer contact with statistics and econometrics in ways that can address practical econometric concerns about model misspecifications and selections of prior probabilities.

Although they proceeded differently than we do, {cite}`Chamberlain:2020`, {cite}`Cerreiaetal:2013`, and {cite}`denti2022model` studied related issues. {cite}`Chamberlain:2020` emphasized that likelihoods and priors are both vulnerable to uncertainties. His ultimate focus was on uncertainty about predictive distributions that are constructed by integrating likelihoods with respect to priors. Our paper instead formulates a decision theory with distinct uncertainties about priors and likelihoods. {cite}`Cerreiaetal:2013` (section 4.2) provide a rationalization of the smooth ambiguity preferences proposed by {cite}`KlibanoffMarinacciMukerji:2005` that includes likelihoods and priors as components. {cite}`denti2022model` extend this approach by using an axiomatic revealed preference approach to deduce an implied parameterization of a likelihood function. But neither of those papers sharply distinguishes prior uncertainty from concerns about possible model misspecifications, something that we want to do. We formulate concerns about model misspecifications as uncertainty about likelihoods.

We assemble concepts and practical ways of modeling risks and concerns about model misspecifications from statistics, robust control theory, economics, and decision theory. We align definitions of statistical models, uncertainty, and ambiguity with ideas from decision theories that build on {cite}`Anscombe_Aumann`'s way of representing subjective and objective uncertainties. We connect our analysis to econometrics and robust control theory by using {cite}`Anscombe_Aumann` states as parameters that index alternative statistical models of random variables that affect outcomes that a decision maker cares about. By modifying {cite}`Gilboaetal:2010`, {cite}`Cerreiaetal:2013`, and {cite}`denti2022model`, we show how to use variational preferences to represent uncertainty about priors and concerns about statistical model misspecifications.

Some "behavioral" models in economics and finance assume expected utility preferences in which an agent's subjective probability differs systematically from probabilities that govern the data.[^behavioral] This literature also contains discussions of differences among agents in their confidence in their view of the world. Lack of confidence can take different forms under different notions of uncertainty. Preference structures that we describe in this paper allow us to formalize different degrees of "confidence" both about details of specifications of particular statistical models and about subjective probabilities to attach to alternative statistical models. Our representations of preferences provide ways to characterize degrees of confidence in terms of perceived statistical plausibilities.[^confidence]

[^deepuncertainties]: Deep uncertainties are defined and discussed by {cite}`hallegatte2012investment`, {cite}`maier2016uncertain`, {cite}`marchau2019decision`, and {cite}`rising2022missing`.

[^econometricians]: Econometricians who explicitly confronted model uncertainty include {cite}`onatski2002robust`, {cite}`brock2003west`, {cite}`stock2006forecasting`, {cite}`brock2007model`, {cite}`del2009monetary`, {cite}`christensen2018dynamic`, {cite}`christensen2019counterfactual`, {cite}`christensen2020robust`, {cite}`andrews2021model`, and {cite}`bonhomme2018minimizing`. {cite}`Chamberlain_2000,Chamberlain_2001` and {cite}`ho:2023` used a post Wald-Savage decision theory of {cite}`GilboaSchmeidler:1989` to confront model uncertainty in his econometric work.

[^likelihoodspriors]: The term likelihood can have multiple meanings. We shall use it to represent a probability density of prize-relevant outcomes, which we refer to as repercussions, conditioned on parameters. Distinguishing likelihood functions from subjective priors is fundamental to Bayesian formulations of statistical learning. See {cite}`deFinetti`, who recommended exchangeability as a more suitable assumption than iid (independent and identically distributed) to model situations in which a decision maker wants to learn. Putting subjective probabilities over parameters that index likelihood functions for iid sequences of random vectors generates exchangeable sequences of random variables.

[^behavioral]: We put "behavioral" in quotes to emphasize that most economic models are about agents' behaviors, including models that impose the rational expectations and common knowledge assumptions that "behavioral" economists want to drop. "Behavioral" economics sometimes means work that is linked more or less informally to psychology.

[^confidence]: Although we provide no formal links to psychology here, we think that a promising research plan would explore connections between so-called behavior distortions and the inferential challenges that economic decision makers confront. As is often assumed in behavioral models, degrees of confidence could differ across economic agents.

**Objects and Interpretations**

Our decision maker knows a parameterized family of probability distributions $\tau(w | \theta)d\upsilon(w),$ where $w \in W$ is a realization of a random vector or "repercussion" that he cares about, $\theta \in \Theta$ is a vector of parameters, and $d\upsilon(w)$ is a measure over $W.$ A realization of $w$ can play two possible roles. It can represent an outcome over which the decision maker has preferences, and it can capture data available to help the decision maker shape decisions. The decision maker has preferences over a set of prize rules, each of which we represent as a function $\gamma: W \times \Theta \rightarrow X$, where $x \in X$ is a "prize" that he cares about. In our featured examples, for parameter vector $\theta \in \Theta$, the prize rule $\gamma(w | \theta)$ determines the decision maker's exposure to an uncertain random vector that has a realization expressible in terms of $w \in W$ A set of $\gamma$'s describes prize rules under consideration. In forecasting problems of a type common in time series statistics and econometrics, the prize can depend directly on the error in forecasting a component of $w$ and the forecast rule depends on another component of $w.$ While forecasting problems are interesting in their own right, in many applications, forecasts are intermediate inputs into outcomes of ultimate interest to the decision maker. Examples that appear in [Preliminaries](sec:prelim) illustrate a range of applications.

````{prf:remark}
:label: rmk:robuststat

We use three components from decision theory, namely, i) states, ii) acts, and iii) prizes, in some new ways. We follow {cite}`Anscombe_Aumann` by defining consequences as lotteries over prizes. An act maps states into consequences. Preferences are defined over acts. In the static setup of this paper, we take the state to be parameters of a statistical model. That distinguishes our formulation from many other applications of {cite}`Anscombe_Aumann`. For example, decision theorists who connect their work to revealed preference theory typically want states that are "verifiable". We are instead interested in a typical econometric situation in which parameters of statistical models remain hidden forever. For us, parameter uncertainty is central, so it is important that parameter vector $theta$ be included as a component of the state.[^paramimportance]

{cite}`Gilboaetal:2010` and {cite}`Cerreiaetal:2013` introduced parameterized models as a family of primitive probabilities that a decision maker cares about. {cite}`Cerreiaetal:2013` in effect consider an expanded state space $(w, theta)$ that includes both repercussions with realization $w$ and parameters $theta$ and then take a *model* to be a conditional distribution over $(W, {\mathfrak W})$ given $theta$.[^dynkinspace] Consistent with the framework of {cite}`Gilboaetal:2010`, {cite}`Cerreiaetal:2013` showed that a family of models induces a partial ordering according to which an act is preferred to another act if it is preferred under all models in the family.

Relative to {cite}`Cerreiaetal:2013` and many other applications of the {cite}`Anscombe_Aumann` framework, we use lotteries in a more essential way. {cite}`Anscombe_Aumann` interpret lotteries as "roulette wheels" with known (objective) probabilities, in contrast to "horse races" with unknown (subjective) probabilities. Many authors used an {cite}`Anscombe_Aumann` setup as a vehicle to extend {cite}`VonneumannMorgenstern:1944` preferences defined over lotteries to more general settings that can include subjective uncertainty. In our formulation, the random vector $W$ induces a probability distribution that according to a particular {cite}`Anscombe_Aumann` act implies a particular lottery that can depend on parameters of a statistical model. We represent a family of models as a family of probability distributions indexed by an unknown parameter vector. Parameter vectors can reside in a finite set, a manifold of possible values, or even an infinite dimensional set. With correct statistical models (i.e., likelihoods), each model induces a "roulette wheel" lottery. The possibility of misspecified likelihoods leads us to want a counterpart to an Anscombe-Aumann lottery with unknown probabilities. Our extension of the {cite}`Anscombe_Aumann` framework lets us distinguish robustness to misspecification of each member of a collection of substantively motivated "structured" statistical models from robustness to the choice of a prior distribution to put over those statistical models. We formulate preferences that express distinct concerns about both types of robustness.

To motivate their axioms, {cite}`MaccheroniMarinacciRustichini:2006b` and {cite}`Strzalecki:2011` used {cite}`HansenSargent:2001`'s stochastic formulation of a robust control problem. We use our {cite}`Anscombe_Aumann` formulation to show that the axioms of {cite}`MaccheroniMarinacciRustichini:2006b` and {cite}`Strzalecki:2011` actually express prior uncertainty and not the model misspecification concerns that had originally motivated {cite}`HansenSargent:2001`. We go on to show how, by using an appropriate ambiguity index or "cost" function, we can use the variational preferences of {cite}`MaccheroniMarinacciRustichini:2006b` to express concerns about robustness both to statistical model misspecification and to prior selection, including priors meant to represent "nonparametric Bayesian" methods.

[^paramimportance]: Stephen Stigler showed us a short working paper by {cite}`Savage:1952` entitled "An Axiomatic Theory of Reasonable Behavior in the Face of Uncertainty," a prolegomenon to the axiomatic structure presented in {cite}`Savage:1954`. {cite}`Savage:1952` wrote this: "The set S represents the conceivable states, or descriptions of the world, or milieu, with which the person is concerned $\ldots$ " We think of parameter values or model selection indicators as presenting a "description of the world."

[^dynkinspace]: {cite}`Cerreiaetal:2013` deploy a "Dynkin space" and an associated sigma algebra of events. Their conditioning on those events is a counterpart to our conditioning on a model. As an alternative, {cite}`denti2022model` used an axiomatic approach to define a parameterized set of models. While both approaches are interesting, we suppose that models can have scientific or other sources from outside the specific decision problem. In this, we follow {cite}`HansenSargent:2022decision` who refer to such models as "structured models."
````
(sec:prelim)=
## Preliminaries

Following {cite}`GilboaSchmeidler:1989` and {cite}`MaccheroniMarinacciRustichini:2006b`, we adopt a version of the framework of
{cite}`Anscombe_Aumann` described by {cite}`Fishburn_70`: $(\Theta, {\mathfrak G})$ is a measurable space of potential *states*, $(X, {\mathfrak X})$ is a measurable space of potential *prizes*, $\Pi$
is a set of probability measures over states, and $\Lambda$ is a set of probability measures over prizes.[^anscombeAumann] For each
$\pi \in \Pi$, $(\Theta, {\mathfrak G}, \pi )$ is a probability space and for each $\lambda \in \Lambda$, $(X, {\mathfrak X}, \lambda)$ is a probability space. Let ${\mathcal X}$ denote an event in ${\mathfrak X}$ and ${\mathcal G}$ denote an event in ${\mathfrak G}$.

[^anscombeAumann]: For a discussion of the Anscombe-Aumann setup, see {cite}`Kreps_notes`, especially chapters 5 and 7.

````{prf:definition}
An **act** is a ${{\mathfrak G}}$ measurable function $f : \Theta \rightarrow \Lambda$.
````

For a given $\theta$, $df(x | \theta)$ denotes integration with respect to the probability measure $f(\theta) \in \Lambda,$ which is a lottery over possible prizes $x \in X$.[^basicsetup] For a **given** probability measure $\pi \in \Pi$, $\int_\Theta df(x | \theta) d\pi(\theta)$ is a two-stage lottery over prizes, with a lottery over states $\theta$ being described by $\pi$ and another lottery over prizes $x \in X$ being described by $df(x | \theta)$, which depends on the outcome $\theta$ from the other lottery. We shall introduce uncertainty about the probability measure $\pi.$

As mentioned in [section](sec:intro), we interpret objects in the {cite}`Anscombe_Aumann` formulation in ways that help us as statisticians/econometricians. We interpret a state $\theta$ as pinning down one among a set $\Theta$ of probability models that a decision maker regards as possible. A decision maker makes a decision (i.e., "chooses an {cite}`Anscombe_Aumann` act") that generates a probability distribution over outcomes that he/she cares about, i.e., over {cite}`Anscombe_Aumann` prizes $x \in X$.

We use {cite}`Anscombe_Aumann` acts to represent alternative conditional distributions of repercussions and prize rules. An action or decision $\delta \in \Delta$, which is distinct from an {cite}`Anscombe_Aumann` act, can be a vector of real numbers or, more generally, a function that is defined on appropriate spaces. A choice of $\delta$ can influence the distribution of repercussions conditioned on the parameter vector. It can also alter the exposure of the prize to repercussions. We represent a decision maker's exposure to repercussions with prize rules $\gamma_\delta$ that are Borel measurable functions that map $W$ into prizes in $X$. We represent the influence of $\delta$ on the distribution of repercussions by a conditional probability measure represented as a density $\tau_\delta( \cdot | \theta)$ with respect to a Borel measure $\upsilon$ on $(W, {\mathcal W}).$ A $\theta \in \Theta$ implies a probability measure

```{math}
\tau_\delta(w | \theta) d\upsilon(w).
```

This formulation is convenient for applied statisticians because for each $\delta,$ a parameterized family $\tau_\delta( \cdot | \theta)$ can define a manifold of likelihoods indexed by a vector of unknown parameters $\theta \in \Theta$. For a given decision $\delta$, $(\gamma_\delta, \tau_\delta)$  
induces a lottery over $X$ conditioned on $\theta,$ and hence an {cite}`Anscombe_Aumann` act that can be represented with $df(x | \theta).$ As we vary decisions $\delta \in \Delta$, we trace out a collection of such acts. A particular decision problem defines both conditional distributions $\tau_\delta d\upsilon$ and prize rules $\gamma_\delta$ for alternative decisions $\delta$. Together, they delineate a collection of {cite}`Anscombe_Aumann` acts.  

[^basicsetup]: We borrow our basic setup from {cite}`Massimo_Simone_2021`. Following the leads of de Finetti and Savage, formulations of max-min expected utility and variational preferences initially worked within a tradition in decision theory under uncertainty that restricted probabilities to be finitely additive. However, much of probability theory routinely imposes countable additivity. It simplifies our presentation.

````{prf:remark}
We can expand the collection of acts by randomizing decisions. Given two decisions $\delta_1$ and $\delta_2,$ a randomized rule chooses decision $\delta_1$ with probability $\alpha$ and $\delta_2$ with probability $1 - \alpha.$ Since each decision induces an {cite}`Anscombe_Aumann` act, the randomized decision is a convex combination of the two induced acts.[^randomization]

[^randomization]: In some special cases, the set of acts induced by decisions may itself be convex. In this case, the randomization of decisions merely replicates the collection of induced acts.
````

We consider several canonical examples.

````{prf:example}
:label: remark:motivation
For some situations, it suffices to let $Delta$ be a Borel set of a finite-dimensional Euclidean space and for the conditional distribution $tau$ not to depend on $delta$. For example,
$delta$ could be a particular portfolio of assets whose random return is exposed to a repercussion vector in a particular way. A choice of a portfolio does not affect the joint distribution of returns on individual assets, but it does influence the return on a portfolio of those component assets.
````

````{prf:example}
:label: example:motivation2
In stochastic optimal control problems like those studied by {cite}`bertsekas1976dynamic`, a decision maker chooses a "control" that affects the distribution of a repercussion, which in this example takes the form of a next-period state vector. For instance, in linear-quadratic Gaussian optimal control problems, often referred to as "linear regulator" problems, this effect shows up in a mean conditioned on a current state. For example, a repercussion vector $w$ obeys:
```{math}
w = {\mathbb A} + {\mathbb B} \delta + {\mathbb C} \epsilon,
``` 
where the probability distribution over $\epsilon$'s is a standard, multivariate normal. The vector ${\mathbb A}$ and the matrices ${\mathbb B}$ and ${\mathbb C}$ depend on parameters in $\Theta$.^[^current_state] Suppose that a controller who chooses $\delta$ knows parameters only up to an uncertain subjective distribution.^[^subjective_distribution] Think of decision $\delta$ as a current period control vector in the sense of {cite}`bertsekas1976dynamic`. The conditional distribution $\tau$ depends on $\delta$: $w$ is distributed as multivariate normal with conditional mean ${\mathbb A} + {\mathbb B} \delta $ and conditional variance ${\mathbb C} {\mathbb C}^\top.$ In typical optimal linear regulator control theory problems, prize rules depend on the vector $(w^\top, \delta^\top)$ with a utility function that is the negative of the quadratic form in this vector, for example,  $- w^\top {\mathbb R} w - \delta^\top {\mathbb Q} \delta$ where $R$ and $Q$ are positive semidefinite matrices. The linear regular problem is an example of a much larger class of stochastic control problems. For simplicity, we formulate it as a static problem.^[^recursive_formulation] The {prf:ref}`remark:motivation` portfolio choice problem is a special case of this stochastic control problem in which repercussions are the returns and a control vector of portfolio weights does not influence repercussions.
````

[^current_state]: The vector ${\mathbb A}$ absorbs the current state. In a standard optimal linear regular problem, the controller knows $({\mathbb A}, {\mathbb B}, {\mathbb C})$.

[^subjective_distribution]: This distribution might depend on past information.

[^recursive_formulation]: More generally, it would be input into a recursive formulation.

````{prf:example}
:label: remark:Ferguson

To build bridges to mathematical statistics, we extend a setup that {cite}`Ferguson:1967` used to analyze learning from data. We again posit a family of densities $tau(w | \theta)$ for a repercussion vector whose realizations are denoted by $w$'s, and where a parameter vector $\theta$ is a "true state of nature". In this example, the decision $\delta$ does not affect $\tau$. We can represent the outcome of what {cite}`Ferguson:1967` calls a statistical experiment as a realization of a random vector $y = \zeta(w)$ that contains information about $w$. This information can be a signal that is correlated with the "prize" of ultimate interest. Let a decision $\delta$ be a measurable function that maps observations $y$ from the statistical experiment into a set of what {cite}`Ferguson:1967` calls *actions*.^[^footnote_decision_rule] In this way {cite}`Ferguson:1967` allows what he calls an "action" -- our decision $\delta$ -- to depend on an observation $y$.

We constrain a prize rule $\gamma_\delta$ to satisfy:^[^footnote_prize_rule]

```{math}
:label: prize_rule
\gamma_\delta (w ) = \Psi[ \delta \circ \zeta (w), w ] .
```

{cite}`Ferguson:1967`'s actions are not {cite}`Anscombe_Aumann` *acts*. To capture {cite}`Ferguson:1967`'s setup, each prize rule $\gamma_\delta$ implies a probability distribution for a prize conditioned on $\theta$ that is induced by $\tau(w | \theta) d\upsilon(w).$ This probability distribution for a prize conditioned on $\theta$ is an {cite}`Anscombe_Aumann` act.

By allowing decisions to depend on data that is observed at an intermediate date, the {cite}`Ferguson:1967` formulation allows a richer collection of possible decisions and nests our Remark {prf:ref}`remark:motivation` formulation as a special case. It can include problems that seek to construct robustly optimal forecasts from historical data. More generally, we are interested in decision problems for which forecasting is an input but not the ultimate goal.

````

Although the problem posed in Example {prf:ref}`remark:Ferguson` is static, it can be reinterpreted as a three-stage or three-period decision problem. A decision rule that is chosen at an initial period zero depends on information about the repercussion that will be revealed in a first stage. The decision maker can condition on this information when choosing his exposure to the repercussion with realization $w$. The repercussion itself is fully realized at the end of stage two. In this formulation, potential likelihood misspecifications affect the decision maker's inferences about the prize distribution in stage one. As posed, this is an *ex-ante* decision problem in which a decision rule, $\delta$, is chosen at period $0$. In contrast, we can view Examples {prf:ref}`remark:motivation` and {prf:ref}`example:motivation2` as *ex post* problems in which the "prior" implicitly conditions on current and past data as does the "decision". As often happens, the timing protocol matters. When a decision maker chooses sequentially, the distinction between priors and posteriors can become obscured when an end of period $j-1$ posterior becomes a period $j$ prior. In a recursive formulation of a dynamic decision problem, concerns about robustness of priors-posteriors can recur in stage-specific components within a multi-stage interpretation of a decision problem. By design, our general formulation invites dynamic extensions.

[^footnote_decision_rule]: For {cite}`Ferguson:1967`, $\delta$ viewed as a function of $y$ is a *decision rule* distinct from our *prize rule*.
[^footnote_prize_rule]: {cite}`Ferguson:1967`'s formulation of the problem introduces a loss function that for us would be the negative of the expectation of a utility function conditioned on $(Y, \theta)$ under the distribution implied by $\tau( w | \theta) d\upsilon(w)$ and $\zeta`.

````{prf:example}
:label: ex:estimation

It is common in econometrics and statistics to pose a decision problem as a parameter estimation problem that supposes that prizes are deviations between an estimator $\delta(y)$ and a function $\chi(\theta)$ of the parameter vector. To capture this, we can let

```{math}
w = \begin{bmatrix} \chi(\theta) \\ y \end{bmatrix}
```

and add a degenerate equation to the $\tau$ dynamics that describes how we construct the first component of $w$. This approach seems to be shorthand for something deeper. Typically, decisions of interest can be expressed in terms of outcomes with probability distributions that depend on the unknown parameter as in the other examples that we mention.
````

In what follows, a decision maker's *prior* over possible statistical models indexed by $\theta$ is a probability measure $\pi \in \Pi$.

````{prf:remark}
The collection of {cite}`Anscombe_Aumann` acts is typically much larger than the set of acts that can be induced by an available pair, $(\gamma_\delta, \tau_\delta)$ for $\delta \in \Delta$ as implied by alternative decisions. We know that the axioms invoked in this paper apply to preferences over the full collection of {cite}`Anscombe_Aumann` acts. While the randomization of decisions described previously enlarges the set of {cite}`Anscombe_Aumann` acts by including the convex hull of the set of acts induced by prize rules, in general that device does not construct the full set of {cite}`Anscombe_Aumann` acts. We recognize that judging the plausibility or "self-evident quality" of the axioms that we impose would require extending the set of the acts to be studied beyond the set of induced by the potential actions within a "substantive decision model" even if allow randomization of the decisions.
````

````{prf:remark}
This is a note of us but will eventually have to be clarified. Alternatively, we might start with a set of acts of the form $f_h$ defined as below that are convex. Thus form $f_h$ and $g_k$ as two such acts. Then $f_{\alpha h + (1-\alpha) k}$ is another act, apparently distinct from the convex combination described previously. What do we say about the preferences between this apparently distinct objects? We will have to say something.
````

Let $\mathcal{A}$ be the set of all acts. Each act $f \in \mathcal{A}$ implies lotteries $f(\theta)$ for each $\theta \in \Theta.$ Two collections of acts will interest us, a set $\mathcal{A}_o$ that lets us represent objective uncertainty and another set $\mathcal{A}_s$ that {cite}`Anscombe_Aumann` used to express subjective uncertainty. Formally, let $\mathcal{A}_o \subset \mathcal{A}$ denote the collection of all *constant* acts where a constant act maps all $\theta \in \Theta$ into a unique lottery over prizes $x \in X$. Constant acts express objective uncertainty because they do not depend on the parameter $\theta$. Absence of dependence means that the probability distribution $\pi \in \Pi$ over states plays no role in shaping an ultimate probability distribution over prizes. A constant act constructed from a prize rule $\gamma$ could emerge as follows. Suppose that some component of $W$ has a known distribution independent of $\theta$ and that $\gamma$ depends only on this component. Such limited dependence implies an act that is independent of $\theta$. The collection $\mathcal{A}_s$ consists of acts, each of which delivers a unique prize for each $\theta$. We let $s(\theta) \in X$ denote an act in $\mathcal{A}_s$.[^acts_footnote] We use a probability distribution $\pi \in \Pi$ over states in conjunction with $\mathcal{A}_s$ to express subjective uncertainty.

[^acts_footnote]: Technically, an act in $\mathcal{A}_s$ is a degenerate Dirac lottery with a mass point at $s(\theta)$ that is assigned probability one.

````{prf:remark}
{cite}`Anscombe_Aumann` distinguished "horse race lotteries," represented by acts in $\mathcal{A}_s$, from "roulette lotteries," represented by acts in $\mathcal{A}_o$. [^kreps_notes] {cite}`Anscombe_Aumann` used roulette lotteries with known probabilities to construct subjective probabilities over horse race events.

[^kreps_notes]: See {cite}`Kreps_notes` for more about the distinction.
````

````{prf:remark}
While {cite}`Savage:1954` did not include "objective" lotteries when he rationalized subjective expected utility, his framework allows flexibility in defining both a state and an act. {cite}`Gilboaetal:2020` exhibit the flexibility of a Savage-style state space with a variety of applications and discuss benefits and challenges that this flexibility brings.[^flexibility] There is also flexibility in constructing an act. {cite}`Carreiaetal:2012` take advantage of this flexibility to produce a preference representation for {cite}`Anscombe_Aumann` acts under {cite}`Savage:1954` axioms augmented with risk independence. This representation coincides with the familiar {cite}`Savage:1954` representation for acts in $mathcal{A}_s$ with unique prizes for each state.[^representation]
````

[^flexibility]: They did not specifically discuss the statistical linkages that we explore here.

[^representation]: More generally, their representation includes an additional curvature adjustment much like the smooth ambiguity model. See Proposition 3 in their appendix.

````{prf:remark}
Though not in ours, in other applications of {cite}`Anscombe_Aumann`, the state is $w$, and uncertainty is about a probability distribution to assign to the space $W$. A lottery conditioned on a state adds additional randomness with known probabilities. Preferences are cast in terms of lotteries conditioned on $w$ as well as uncertainty over $W$.
````

{cite}`Anscombe_Aumann` wrote:

> "... anyone who wishes to avoid a concept of physical chance distinct from probability may reinterpret our construction as a method of defining more difficult probabilities in terms of easier ones."

This environment is often interpreted as containing both "subjective" and "objective" uncertainties (for example, see {cite}`Strzalecki`).[^AnscombeAumannFootnote]

[^AnscombeAumannFootnote]: {cite}`Anscombe_Aumann` distinguished "horse lotteries" with unknown probability distributions from "roulette lotteries" having known probability distributions. See {cite}`Kreps_notes` for more about the distinction. Thus, because the (unknown) "prior" distribution over states $\theta$ belongs to a nontrivial set $\Pi$, there is "subjective" uncertainty about $\theta$. But given $\theta$, $f(\theta)$ is a known distribution over prizes $x$, so there is "objective" uncertainty about $x \in X$. According to these distinctions, we can imagine:
-  Purely objective uncertainty that takes the form of a constant act, an act that delivers the **same** lottery $\lambda \in \Lambda$ for all $\theta \in \Theta$.
-  Purely subjective uncertainty that takes the form of an act, $f$, that for each $\theta$, for which $f(dx | \theta)$ is a degenerate Dirac lottery that pays prize $h(\theta) \in X$ with probability one and the prize  $h(\theta)$  can depend on $\theta$.

We shall often construct a new act from initial acts $f$ and $g$ by using a probability $\alpha \in (0,1)$ to form a mixture

```{math}
\left[ \alpha f + (1 - \alpha) g \right] (\theta) = \alpha f(\theta) + (1 - \alpha) g(\theta) \in \Lambda \hspace{.5cm} \forall \theta \in \Theta .
```

We shall use instances of our {cite}`Anscombe_Aumann` framework to describe a) a Bayesian decision maker with a unique prior over a set $\Theta$ of statistical models, b) a decision maker who knows a set $\Theta$ of statistical models and who copes with **ambiguity** about those models by considering prospective outcomes under a set of priors $\Pi$ over those statistical models, c) a decision maker with concerns that a single known statistical model $\theta$ is **misspecified** by using a statistical discrepancy measure to delineate unknown models surrounding that known model, and d) a decision maker with ambiguity and concerns about model misspecifications.

### Preferences

To represent a decision maker's preferences over acts, we use $\sim$ to mean indifference, $\succsim$ a weak preference, and $\succ$ a strict preference. Throughout, we assume that preferences are non-degenerate (there is a strict ranking between two acts), complete (we can compare any pair of acts), and transitive ($f \succsim g$ and $g \succsim h$ imply $f \succsim h$). We also impose an Archimedean axiom that provides a form of continuity.[^footnoteArch]

[^footnoteArch]: The Archimedean axiom states: let $f,g,h$ be acts in ${\mathcal A}$ with $f \succ g \succ h$. Then there are $0 < \alpha_1 < 1$ and $0 < \alpha_2 < 1$ such that $\alpha_1 f + (1-\alpha_1) h \succ g \succ \alpha_2 f + (1-\alpha_2) h$. See {cite}`Herstein_Milnor` for an alternative formulation of a continuity axiom.

A **finite signed measure** on the measurable space $(X, \mathfrak{X})$ is a finite linear combination of probability measures that resides in a linear space ${\widehat \Lambda}$ that contains $\Lambda$.

### Objective probability

By analyzing preferences over the constant acts ${\mathcal A}_o$, we temporarily put aside attitudes about ambiguity and model misspecification and focus on objective uncertainty (sometimes called "risk"). There is a unique probability $\lambda \in \Lambda$ associated with every act $f \in {\mathcal A}_o$ and a unique act in $ {\mathcal A}_o$ associated with every $\lambda \in \Lambda$. We define a restriction $ \succ_\Lambda$ of the preference order $ \succ$ to the space of constant acts $f \in {\mathcal A}_o$ by

```{math}
\lambda \succ_\Lambda  \kappa \iff f \succ g
```

where $\lambda$ is the probability generated by act $f \in {\mathcal A}_o$ and $\kappa$ is the probability distribution generated by act $g \in {\mathcal A}_o$.

To represent preferences $\succ_\Lambda$, we follow {cite}`VonneumannMorgenstern:1944`, who imposed the following restriction:[^footnoteNonDeg]

[^footnoteNonDeg]: Completeness, transitivity, and the Archimedean axiom carry over directly from $\succ$ to $\succ_\Lambda$, but not necessarily non-degeneracy. Our presentation below presumes non-degeneracy of $\succ_\Lambda`.


````{prf:axiom} 
:label: ax:ind0
(Independence) If $f, g, h \in {\mathcal A}_o$ and $\alpha \in (0,1)$, then

```{math}
f \succsim g \Rightarrow \alpha f + (1-\alpha) h \succsim  \alpha g + (1 - \alpha) h .
```
````

The approach delivers an expected utility representation of preferences over constant acts: there exists a utility function $u: X \rightarrow \mathbb{R}$ such that for $f, g \in {\mathcal A}_o$

```{math}
f \succsim g \iff  U(f) \ge U(g)
```

where

```{math}
U(f) = \int_X u(x) d \lambda(x)
```

and $\lambda \in \Lambda$ is the probability distribution generated by constant act $f$.

### Subjective probability

To construct subjective expected utility preferences, we extend an expected utility representation of $\succ_\Lambda$ to a representation of preferences $\succ$ on the set ${\mathcal A}$ of all acts. We impose restrictions on $\succ$ in the form of two axioms.

The first extends the independence axiom to the set of all acts:


````{prf:axiom} 
:label: ax:ind
(Independence) If $f, g, h \in {\mathcal A}$ and $\alpha \in (0,1)$, then:

```{math}
f \succsim g \Rightarrow \alpha f + (1-\alpha) h \succsim  \alpha g + (1 - \alpha) h .
```
````

The second is:


````{prf:axiom} 
:label: ax:postive
(Monotonicity) For any $f, g \in {\mathcal A}$ such that $f(\theta) \succsim_\Lambda g(\theta)$ for each $\theta in \Theta$, $f \succsim g$.
````

This approach leads to the following representation of preferences over acts $f \in {\mathcal A}$

```{math}
f \succsim g \iff \int_{\Theta} \left[\int_X u(x) d f(x | \theta) \right] d \pi(\theta) \ge \int_\Theta \left[ \int_X u(x) d g(x | \theta) \right] d \pi(\theta),
```

where the probability measure $\pi$ describes subjective probabilities.

### Max-min Expected Utility

To construct a decision maker who has max-min expected utility preferences, {cite}`GilboaSchmeidler:1989` replaced Axiom {prf:ref}`ax:ind` with the following two axioms:


````{prf:axiom} 
:label: ax:gs1
(Certainty Independence) If $f, g \in {\mathcal A}$, $h \in {\mathcal A}_o$, and $\alpha \in (0,1)$, then:

```{math}
f \succsim g \iff \alpha f + (1-\alpha) h \succsim  \alpha g + (1 - \alpha) h.
```
````


````{prf:axiom} 
:label: ax:gs2
(Uncertainty Aversion) If $f, g \in {\mathcal A}$ and $\alpha \in (0,1)$, then:

```{math}
f \sim g \Rightarrow  \alpha f + (1-\alpha)g \succsim f.
```
````

````{prf:example}
Suppose that $\Theta = \{\theta_1, \theta_2\}$ and consider lotteries $\lambda_1$ and $\lambda_2$. Let act $f$ be lottery $\lambda_1$ if $\theta = \theta_1$ and lottery $\lambda_2$ if $\theta = \theta_2$. Let act $g$ be lottery $\lambda_2$ if $\theta = \theta_1$ and lottery $\lambda_1$ if $\theta = \theta_2$. Suppose that $f \sim g$. Axiom {prf:ref}`ax:gs2` allows a preference for mixing the two acts. If, for instance, $\alpha = \frac{1}{2}$, the mixture is a constant act with a lottery $\frac{1}{2} \lambda_1 + \frac{1}{2} \lambda_2$ that is independent of $\theta$. We think of mixing as reducing the exposure to $\theta$ uncertainty. In the extreme case, setting $\alpha = \frac{1}{2}$, for example, completely eliminates effects of exposure to $\theta$ uncertainty.
````

By replacing Axiom {prf:ref}`ax:ind` with Axioms {prf:ref}`ax:gs1` and {prf:ref}`ax:gs2`, {cite}`GilboaSchmeidler:1989` obtained preferences described by
```{math}
:eq:GS101
f \succsim g \iff \min_{\pi \in \Pi_c} \int_\Theta \left[\int_X u(x) df(x \mid \theta)\right] d\pi(\theta) \ge \min_{\pi \in \Pi_c} \int_\Theta \left[\int_X u(x) dg(x \mid \theta)\right] d\pi(\theta)
```
for a convex set $\Pi_c \subset \Pi$ of probability measures. An act $f(\theta)$ is still a lottery over prizes $x \in X$ and, as in representation \eqref{representation}, for each $\theta$, $\int_X u(x) df(x \mid \theta)$ is an expected utility over prizes $x$.
Evidently, expected utility preferences \eqref{representationfishburn} are a special case of max-min expected utility preferences {eq}`eq:GS101` in which $\Pi_c$ is a set with a single member.

(sec:varpref1)=
## Variational preferences

{cite}`MaccheroniMarinacciRustichini:2006b` relaxed certainty independence Axiom {prf:ref}`ax:gs1` of {cite}`GilboaSchmeidler:1989` to obtain preferences with a yet more general representation that they called variational preferences.


````{prf:axiom}
:label: ax:mmr2
(Weak Certainty Independence) If $f, g \in \mathcal{A}$, $h, k \in \mathcal{A}_o$, and $\alpha \in (0,1)$, then
```{math}
\alpha f +(1-\alpha) h \succsim \alpha g +(1-\alpha) h \Rightarrow \alpha f +(1-\alpha) k \succsim \alpha g +(1-\alpha) k
```
````

Axiom \ref{ax:mmr2} considers only acts that are mixtures of constant acts that can be represented with a single lottery. The axiom states that altering the constant act from $h$ to $k$ does not reverse the decision maker's preferences. The same $\alpha$ appears in all three acts being compared. This axiom imparts to preferences a smooth tradeoff between separate contributions that come from an expected utility, on the one hand, and from statistical uncertainty, on the other.

{cite}`MaccheroniMarinacciRustichini:2006b` showed that preferences that satisfy the weaker Axiom \ref{ax:mmr2} instead of Axiom \ref{ax:gs1} are described by
```{math}
:eq:MMR101
f \succsim g \iff \min_{\pi \in \Pi} \int_\Theta \left[\int_X u(x) df(x \mid \theta)\right] d\pi(\theta) + c(\pi) \ge \min_{\pi \in \Pi} \int_\Theta \left[\int_X u(x) dg(x \mid \theta)\right] d\pi(\theta) + c(\pi)
```
where $u$ is uniquely determined up to a linear translation and $c$ is a convex function that satisfies $\inf_{\pi \in \Pi} c(\pi) = 0$. Smaller convex functions, $c$, express more aversion to uncertainty. The convex function $c$ in variational preferences representation \eqref{eq:MMR101} replaces the restricted set of probabilities $\Pi_c$ that appears in the max-min expected utility representation {eq}`eq:GS101`.

(sec:statdivergasc)=
## Scaled statistical divergences as $c$ functions

*Scaled statistical divergences* give rise to convex $c$ functions that especially interest us. We use such divergences in two ways, one for distributions over $(W, {\mathfrak W})$, another for distributions over $(\Pi, {\mathfrak G} )$. We construct statistical divergences for these two situations in similar ways.

We first consider repercussion distributions over $(W, {\mathfrak W})$. Consider a family of probabilities represented as densities with respect to $\upsilon$:

```{math}
:label: density_set
{\mathcal L} := \left\{ \ell \ge 0 : \int   \ell(w) d\upsilon(w) = 1 \right\}
```

For a baseline density $\ell_o$, a *statistical divergence* is a convex function $D(\ell \mid \ell_o)$ of probability measures $\ell(w) d \upsilon(w)$ that satisfies

- $D(\ell \mid  \ell_o) \geq 0$
- $D(\ell \mid  \ell_o) = 0$ implies $\ell = \ell_o$

Given $\ell_o,$ write:

```{math}
\ell (w) d\upsilon(w) = m(w) \ell_o(w) d \upsilon(w)  
```

for $m = \frac {\ell}{ \ell_o}$, where we assume $m$ is not infinite with positive $\upsilon$ measure so that the probability measure $\ell$ is absolutely continuous with respect to $\ell_o(w) d\upsilon(w)$.[^divInfiniteM] The set of such densities is convex as is the set of implied relative densities $m$. To define a scaled statistical divergence, we set

```{math}
D(\ell \mid  \ell_o) = \xi \int_W \phi[m(w) ] \ell_o(w) d \upsilon(w),
```

where $\xi > 0,$ and $\phi$ is a convex function defined over the nonnegative real numbers for which $\phi(1) = 0$ and impose $\phi''(1) = 1$ as a normalization. Examples of such $\phi$ functions and the divergences that they lead to are

```{math}
\begin{align*}
\phi(m) &= - \log (m) & \text{Burg entropy}  \\
\phi(m) & = -4  \left( \sqrt{m} - 1 \right) & \text{ Hellinger distance} \\
\phi(m) & = m \log (m) & \text{ relative entropy}\\
\phi(m) & = {\frac 1 2} \left( m^2 - m \right) & \text{quadratic}.
\end{align*}
```

When $\xi = 1$, the divergence, $D,$ is often called a $\phi$ or $f$-divergence. When $\phi(m)  = m \log (m)$ and $\xi = 1$, we obtain relative entropy

```{math}
D_{KL}(\ell | \ell_o) = \int_W m(w)\log[m(w)]  d \upsilon(w)    .
```

Relative entropy is commonly referred to as Kullback-Leibler divergence.

[^divInfiniteM]: For $\ell$'s for which the implied $m$ is infinite with positive $\upsilon$ measure, we define the divergence to be infinity.

````{prf:remark}
Other families of divergences can be used in conjunction with preference representations that follow, for instance, from Bregman and Wasserstein divergences.  
The family $\phi$ or $f$ divergences featured here has very nice duality properties.   
as we will see, duality allows us to make formal connections to the extensive literature on smooth ambiguity.  
Furthermore, these divergences are invariant to one-to-one transformations of the space over which the probability distributions are defined.  In addition,  some members of this family have useful links  to statistical discrimination procedures. The link to likelihood-based statistical discrimination enables statistical constructions that can help us calibrate concerns about robustness.
````

(sec:ourformulation)=
## Basic formulation

We associate a probability measure $\tau(w | \theta)d\upsilon(w)$ parametrized by $\theta \in \Theta$ with a random vector having possible realizations $w$ in the measurable space $(W, \mathfrak{W})$. Consider alternative real valued, Borel measurable functions $\gamma \in \Psi$ that map $w \in W$ into an $x \in X$. Think of $\gamma$ as a prize rule and $\gamma(w)$ as an uncertain scalar prize. For each prize rule $\gamma$, let $d \lambda(x \mid \theta)$ be the distribution of the prize $x$ that is induced by distribution $\tau(w | \theta) d\upsilon(w)$ and the prize rule $\gamma$. The distribution of the prize thus depends both on the prize rule $\gamma(w)$ and the distribution $\tau(w \mid \theta) d\upsilon(w).$ Within this setting, a decision $\delta$ gives rise to a specific pair $(\gamma_\delta, \tau_\delta)$. To avoid cluttering our notation, we will drop the explicit dependence of $\gamma$ and $\tau$ on $\delta$ in much of the following discussion.

(sec:nonknowingprior)=
### Not knowing a prior

Like the robust Bayesian decision maker of {cite}`Berger:1984`, {cite}`Gilboaetal:2010` and {cite}`Cerreiaetal:2013`, our decision maker has multiple prior distributions because he does not trust the baseline prior.[^robustBayesian] We label such distrust of a single prior "model ambiguity." Here we describe a static version of what {cite}`HansenSargent:2020,HansenSargent:2022decision` call structured uncertainty. "Structured" refers to a particular way that we reduce the dimension of a set of alternative models relative to the much larger set considered when we explore likelihood or model misspecification.

A baseline $\pi_o$ anchors a set of priors $\pi$ over which a decision maker wishes to be robust. We describe the set of priors by
```{math}
\pi(d\theta) = n(\theta)\pi_o(d\theta) ,
```
where $n$ is in the set $\mathcal{N}$ defined by:
```{math}
:label: Nset
\mathcal{N} \doteq \left\{ n \ge 0  : n(\theta) \ge 0, \int_\Theta n(\theta)  d\pi_o(\theta) = 1 \right\}.
```
This specification includes a form of "structured" uncertainty in which all models have the same parametric "structure" but in which each is associated with a different vector of parameter values.[^structuredUncertainty] The decision maker is certain about each of the specific models but is uncertain about a prior to put over them.

(sec:notknowingprior1)=
#### Not knowing a prior, I

To express a form of ambiguity aversion, the decision maker uses scaled statistical divergence
```{math}
c(\pi) = \xi \int_\Theta \phi \left[ n(\theta) \right]    d\pi_o(\theta) 
```
and has variational preferences ordered by[^CerreiaetalVariational]
```{math}
:label: eqn:minproblem106
\min_{n \in \mathcal{N} } \int_\Theta \left(  \int_W u[\gamma(w) ]  \tau( w \mid \theta) d\upsilon(w)  \right) n(\theta) d\pi_o(\theta)    + \xi \int_\Theta \phi[n(\theta)] d \pi_o(\theta)   .
```

[^robustBayesian]: See {cite}`Berger:1984` for a robust Bayesian perspective. By applying a {cite}`GilboaSchmeidler:1989` representation of ambiguity aversion to a decision maker who has multiple predictive distributions, {cite}`Cerreiaetal:2013` forge a link between ambiguity aversion as studied in decision theory and the robust approach to statistics.
[^structuredUncertainty]: See {cite}`HansenSargent:2022decision`.
[^CerreiaetalVariational]: See Theorem 4 of {cite}`Cerreiaetal:2013` for their counterpart to this representation.

````{prf:remark}
It is convenient to solve the minimization problem {eq}`eqn:minproblem106` by using duality properties of convex functions. Because the objective is separable in $\theta$, we first compute
```{math}
:label: eqn:legendre102
\phi^*({\sf u} \mid \xi ) = \min_{{\sf n} \ge 0} {\sf u} {\sf n} + \xi \phi({\sf n})
```
where ${\sf u} = \int u[\gamma(w)]\tau(w | \theta)d\upsilon(w) + \eta$, ${\sf n}$ is a nonnegative number, and $\eta$ is a nonnegative real-valued Lagrange multiplier attached to the constraint
$\int_\Theta n(\theta) d \pi_o(\theta) = 1$; $\phi^*({\sf u} \mid \xi )$ is a concave function of ${\sf u}.$[^legendretransform] The minimizing value of ${\sf n}$ satisfies
```{math}
{\sf n}^* = \phi'^{-1} \left(- {\frac {\sf u} \xi} \right) .
```
The dual to the minimization problem on the right side
of {eq}`eqn:minproblem106` is
```{math}
:label: eqn:dualproblem106
\max_\eta \int_\Theta \phi^* \left(u[\gamma(w)] \tau(w \mid \theta) d\upsilon(w) +\eta \right) d \pi_o(\theta) - \eta .
```

[^legendretransform]: The function $- \phi^*(-{\sf u} \mid \xi)$ is the Legendre transform of $\xi \phi({\sf n})$.
````

````{prf:remark}
(Smooth ambiguity preferences)

When statistical divergence is scaled relative entropy, preferences over $\gamma(w)$ are ordered by

```{math}
:label: eqn:risksens_1015
- \xi \log \left[  \int_\Theta \exp \left( - {\frac { \int_W u[ \gamma(w) ]  \tau( w  \mid \theta) d\upsilon(w)} \xi}\right) d \pi_o(\theta) \right],
```

a static version of preferences that {cite}`HansenSargent:2007` used to frame a robust dynamic filtering and control problem.

These preferences are also a special case of the smooth ambiguity preferences that {cite}`KlibanoffMarinacciMukerji:2005` justified with a set of axioms different from the ones we have used here. Furthermore, {cite}`MaccheroniMarinacciRustichini:2006b` and {cite}`Strzalecki:2011` use this formulation to express concerns about model misspecification.[^footnoteStrzalecki]

In contradistinction, the robustness concerns being represented in this subsection are about a baseline prior over known models and not about possible misspecifications of those models.

[^footnoteStrzalecki]: {cite}`Strzalecki:2011` showed that when Savage's Sure Thing Principle augments axioms imposed by {cite}`MaccheroniMarinacciRustichini:2006b`, the cost functions capable of representing variational preferences are proportional to scalar multiples of entropy divergence relative to a unique baseline prior. The Sure Thing Principle also plays a significant role in {cite}`denti2022model`'s axiomatic construction of a parameterized likelihood to be used in {cite}`KlibanoffMarinacciMukerji:2005` preferences.
````

````{prf:remark}
In subsection [Not knowing a model](sec:notknowingmodel), we considered a setting in which the state or "parameter vector" $\theta = m$ and

```{math}
:label: eqn:stateparameter
d \tau( w \mid \theta ) = m(w) d \tau_o(w),
```

where $m \in {\mathcal M}$ and $d \tau_o (w)$ is a baseline distribution. Suppose that instead of proceeding as we did in subsection [Not knowing a model](sec:notknowingmodel) we were to posit a baseline prior $\pi_o$ over a possibly infinite dimensional parameter space ${\mathcal M}$. If we were then to explore prior robustness by considering alternative $n$'s in ${\mathcal N}$ using statistical divergence [section](sec:statdivergasc), we would have to specify $\pi_o$ to be sufficiently "informative" that it would limit the range of alternative distributions $d \tau(w)$ that the decision maker entertains relative to those in the [Not knowing a model](sec:notknowingmodel) formulation. The distinct ways in which the [Not knowing a model](sec:notknowingmodel) and [Not knowing a prior](sec:nonknowingprior) formulations use statistical discrepancies lead to substantial differences in the resulting variational preferences, namely, representation {eq}`eqn:minproblem105` for the [Not knowing a model](sec:notknowingmodel) setting of not knowing the distribution of $d \tau(w)$ and {eq}`eqn:minproblem106` for the [Not knowing a prior](sec:nonknowingprior) formulation of not knowing a prior over a known set of known models $d \tau(w \mid \theta)$. See remark {prf:ref}`rem:sims` for more about this issue.
````

(sec:nonknowingprior2)=
### Not knowing a prior, II

We modify preferences by using a statistical divergence to constrain a set of prior probabilities. The resulting preferences satisfy axioms of {cite}`GilboaSchmeidler:1989`. Consider:

```{math}
:label: Pidef
\Pi = \{ \pi : d\pi(\theta) = n(\theta) d\pi_o(\theta), n \in {\mathcal N}, \int_\Theta \phi[n(\theta)] d \pi_o(\theta) \le \kappa  \}
```
where $\kappa > 0$ pins down the size of the set of priors. Preferences over $\gamma(w)$ are ordered by

```{math}
:label: eqn:minproblem107
\min_{\pi \in \Pi } \int_\Theta \left( \int_W u[ \gamma( w )] \tau( w \mid \theta) d\nu(w) \right) d \pi(\theta).
```

````{prf:remark}
It is convenient to solve the minimization problem on the right side of {eq}`eqn:minproblem107` by using duality properties of convex functions. The minimized objective for problem {eq}`eqn:minproblem107` can again be evaluated using convex duality theory. We now explicitly note the dependence of $\phi^*$ on $\xi$ and write the dual problem as:

```{math}
\max_{\eta, \xi \ge 0}   \int_\Theta \phi^* \left[  \int_W  u[\gamma(w)] \tau(w \mid \theta) d\nu(w)   + \eta  \mid \xi \right] d \pi_o(\theta ) - \eta - \xi \kappa .
```

Maximization over $\xi \ge 0 $ enforces a constraint on the set of admissible priors.
````

````{prf:remark}
Within a setting like that of {prf:ref}`ex:estimation`, {cite}`ho:2023` used another approach to compute robust adjustments to posterior expectations. He used this approach to assess the prior sensitivity of empirical measurements of targets of interest to an investigator. {cite}`ho:2023`'s framework could also be used to define robust preferences defined in terms of posterior expectations. For instance, measurements of interest could be depicted as the maximizer of the negative of an expected loss function of a type common in statistics and econometrics. More formally, {cite}`ho:2023` used relative entropy divergence to restrict a set of priors. He computed expectations conditioned on a signal and minimized over possible implied posterior distributions given a relative entropy constraint over the priors. The minimizing "prior" from this approach typically depends on the signal,[^hoformula] unlike the outcome from solving the *ex-ante problem* described in {prf:ref}`remark:Ferguson`. Dependence of a minimizing prior on the signal like that in {cite}`ho:2023`'s formulation also emerges in some recursive formulations of dynamic problems, a situation that can lead to statistically inadmissible decisions.[^dynamicproblem]
````

[^hoformula]: See formula (2.7) in {cite}`ho:2023`.

[^dynamicproblem]: See {cite}`Epstein_Schneider_2003` for a formulation of a dynamic choice problem under ambiguity aversion that deploys multiple priors recursively. {cite}`HansenSargent:2022decision` described a possible tension between admissibility and dynamic consistency.


(sec:notknowingmodel)=
## Not knowing a likelihood

Instead of being about a prior as the previous subsection, we now suppose that the decision maker's uncertainty is about a likelihood function. We start by supposing that there is a single model that the decision maker fears is misspecified. We then extend the analysis by introducing a parameterized family of probability models that a decision maker thinks might be misspecified.

(sec:notknowingmodel1)=
### A misspecified model

Consider first a single model that might be misspecified. We study a decision maker who knows a parameter $theta_o$. We also fix a decision $delta$, a determinant of $tau$ that we continue to leave implicit in our notation. The decision maker entertains the possible misspecification of

```{math}
tau_o(w) := tau(\cdot \mid \theta_o)
```

in ways that the decision maker cannot precisely describe. But he can say that the alternative models that he is most worried about are statistically close to his baseline model. The presence of too many statistically nearby models would prevent a Bayesian from deploying a proper prior over them. (Later we will compare our approach here to a robust Bayesian approach that requires a family of priors that are mutually absolutely continuous.)

Notice that $\tau_o in {\mathcal L} $ where ${\mathcal L}$ is given by 

```{math}
:label: misspecified_models
m(w) \tau_o(w) 
```

where

```{math}
m(w) := \frac{\ell} {\tau_o}
```

$\ell in {\mathcal L}$ for ${\mathcal L}$ given by {eq}`density_set`. We represent the decision maker's ignorance of specific alternative models by assuming that he entertains a potentially infinite dimensional space ${\mathcal L}$ of what we will call "unstructured" models. A decision maker's expected utility under alternative model $\ell \tau_o(w) d\upsilon(w)$ is

```{math}
:label: eqn:exputilmonly2
\int_W u[\gamma(w)] m(w) \tau_o(w) d\upsilon(w) =  \int_W u[\gamma(w)] \ell (w)  d\upsilon(w).
```

Notice that {eq}`eqn:exputilmonly2` evaluates expected utility for a single choice for $m$.

````{prf:remark}
:label: rem:sims
Specifying a prior over the infinite dimensional space $\mathcal{M}$ brings challenges associated with all nonparametric methods, including "nonparametric Bayesian" methods. A Bayesian prior on an infinite dimensional space such as $\mathcal{M}$ must be more "informative" than is required in finite-dimensional estimation problems.[^sims2010] A related "informativeness" requirement carries over to families of priors that are absolutely continuous relative to a baseline prior. The decision maker in this subsection does not want to entertain priors that are "too informative." In subsection [nonknowingprior1](sec:nonknowingprior1), we describe a decision maker who is concerned about a set of models that is small enough to proceed with a "robust Bayesian" approach with priors over those models that are not "too informative."

[^sims2010]: {cite}`sims2010understanding` critically surveys an extensive statistical literature on this issue. Foundational papers are {cite}`Freedman_1963`, {cite}`SimsAnnals71`, and {cite}`Diaconis_Freedman_1986`.

````

To complete a description of preferences, we require a scaled statistical divergence. We consider alternative probabilities parameterized by entries in $\mathcal{L}$. Under this perspective, a probability model corresponds to a choice of $m \in \mathcal{M}$. Form a scaled divergence measure:

```{math}
:label: eqn:scaledKL99
c( m) = \xi \int_{W} \phi [m(w) ] \tau_o(w) d\upsilon(w)
```
where $\xi > 0$ is a real number. Variational preferences that use {eq}`eqn:scaledKL99` as scaled statistical divergence are ordered by

```{math}
:label: eqn:minproblem105
\min_{m = \frac{\ell}{\tau_o}, \ell \in \mathcal{L}} \left( \int_{W} u[\gamma(w)] m(w) \tau_o(w) d\upsilon(w) + \xi \int_{W} \phi[m(w)] \tau_o(w) d\upsilon(w) \right).
```
This formulation lets a decision maker evaluate alternative prize rules $\gamma(w)$ while guarding against a concern that his baseline model $\tau_o$ is misspecified without having in mind specific alternative models $\tau$. Key ingredients are the single baseline probability $\tau_o$ and a statistical divergence over probability distributions $m(w) \tau_o(w) d\upsilon(w)$.

````{prf:remark}
:label: rem:dual1
As was the case for robust prior analysis, it is again convenient to solve the minimization problem on the right side of {eq}`eqn:minproblem105` by using duality properties of convex functions. Because the objective is separable in $w$, we can first compute

```{math}
:label: eqn:legendre102
\phi^*(\mathsf{u} \mid \xi ) =  \min_{\mathsf{m} \ge 0}  \mathsf{u} \mathsf{m}  + \xi \phi(\mathsf{m})
```
where $\mathsf{u} =u[\gamma(w)] + \eta$, $\mathsf{m}$ is a nonnegative number, and $\eta$ is a nonnegative real-valued Lagrange multiplier that we attach to the constraint 
$\int m(w) \tau_o(w )d\upsilon(w) = 1$; $\phi^*(\mathsf{u} \mid \xi )$ is a concave function of $\mathsf{u}$. The minimizing value of $\mathsf{m}$ now satisfies

```{math}
\mathsf{m}^*  =  \phi'^{-1}  \left(- {\frac {\mathsf{u}}  \xi} \right) .
```
The dual problem to the minimization problem on the right side
of {eq}`eqn:minproblem105` is

```{math}
:label: eqn:dualproblem105
\max_\eta \int_W  \phi^* (u[\gamma(w)] +\eta  ) \tau_o(w)d\upsilon(w)  - \eta .
```
````

````{prf:remark}
:label: rem:dual2
We posed a minimum problem {eq}`eqn:minproblem105` in terms of a set of probability measures on the measurable space $(W, \mathfrak{W})$ with baseline probability $\tau_o(w)d\upsilon(w)$. Since the integrand in the dual problem {eq}`eqn:dualproblem105` depends on $w$ only through the control law $\gamma$, we could instead have used the same convex function $\phi$ to pose a minimization in terms of a set of probability distributions $d\lambda(x)$ with the baseline being the probability distribution over prizes induced $x = \gamma(w)$ with distribution $d\lambda_o(x)$. Doing that would lead to equivalent outcomes. Representations in sections [Preliminaries](sec:prelim) and [Variational Preferences](sec:varpref1) are all cast in terms of induced distributions over prizes. Because control problems entail searching over alternative $\gamma$'s, it is more convenient to formulate them in terms of a baseline model $\tau_o(w)d\upsilon(w)$, as we originally did in subsection [Not Knowing the Model](sec:notknowingmodel).
````

````{prf:remark}
:label: rem:dual3

If we use relative entropy as a statistical divergence, then

```{math}
\phi^*({\sf u} \mid \xi ) =  - \xi  \exp \left( - {\frac {\sf u + \eta } \xi} - 1  \right)
```

and dual problem {eq}`eqn:dualproblem105` becomes[^DupuisEllis]

```{math}
:label: eqn:risksens1002
\max_\eta  - \xi  \int_W  \exp \left[ - {\frac { u[\gamma(w)]  + \eta } \xi} - 1 \right] \tau_o(w) d \upsilon(w)  - \eta = - \xi \log \left(  \int_{W}  \exp \left[ - {\frac { u[\gamma(w)]} \xi}\right]  \tau_o(w) d\upsilon(w) \right) .
```

The minimizing $m$ in problem {eq}`eqn:minproblem105` is

```{math}
:label: eq:n*
m^*(w) =  \frac {\exp\left[ - {\frac {u[\gamma(w)]} \xi}  \right] } { \int_W \exp\left[ - {\frac {u[\gamma(w)] }  \xi} \right] \tau_o(w) d\upsilon(w)}.
```

The worst-case likelihood ratio $m^*$ exponentially tilts a lottery toward low-utility outcomes. {cite}`Bucklew_2004` calls this adverse tilting a statistical version of Murphy's law:

> "The probability of anything happening is in inverse proportion to its desirability."

[^DupuisEllis]: See {cite}`Dupuis_Ellis` for a closely related connection between relative entropy and a variational formula that occurs in large deviation theory.
````

````{prf:remark} 
:label: rem:risksens100 
(Risk-sensitive preferences)

The right side of equation {eq}`eqn:risksens103`, namely,

```{math}
:label: eqn:risksens103
- \xi \log \left[  \int_W \exp \left( - {\frac { u[\gamma(w)]} \xi}\right)  \tau_o(w) d\upsilon(w) \right] ,
```

defines what are known as "risk-sensitive" preferences over control laws $\gamma$. Since a logarithm is a monotone function, these are evidently equivalent to {cite}`VonneumannMorgenstern:1944` expected utility preferences with utility function

```{math}
- \exp\left[ - {\frac {u(\cdot)} \xi} \right]
```

in conjunction with the baseline distribution $\tau_o$ over repercussions. Risk-sensitive preferences are widely used in robust control theory (for example, see {cite}`Jacobson:1973`, {cite}`whittle1990,whittle1996optimal`, and {cite}`PetersenJamesDupuis:2000`).

````

````{prf:remark}
Although our notation suppressed it, the $m$'s in the minimization problem can depend on the decision $\delta$, as dependence that carries over to implied densities, $\ell.$
````

````{prf:example}
:label: ex:discrete

We could say that {eq}`misspecified_models` gives a parameterization of alternative models expressed in terms of $m$ or $\ell$. But since the divergence {eq}`eqn:minproblem105` is not expressed as a divergence in terms of priors over the parameter space, we then could not view preferences {eq}`eqn:minproblem105` as a special case of the robust Bayesian decision theory described in [section](sec:notknowingprior1). With potential misspecifications present, we have deliberately avoided imposing a baseline prior over $\mathcal{L}$.[^divpref]
Instead, each $m$ induces an alternative {cite}`Anscombe_Aumann` roulette wheel. For reasons that will become clear in the next subsection, we think of $m$s as ways to introduce ambiguities about lotteries, disarming the "roulette wheel" analogy.

[^divpref]: Divergence preferences typically are expressed in terms of probabilities distributions over the collection states, which in the present case would be over the set $\mathcal{L}$'s.
````

````{prf:remark}
Preferences that use a relative entropy divergence to capture concerns about model misspecification are often referred to as "multiplier preferences." Because of the different ways that we apply the language of decision theory, the preceding construction of multiplier preferences differs from constructions provided by {cite}`MaccheroniMarinacciRustichini:2006b` and {cite}`Strzalecki:2011`. Specifically, {cite}`MaccheroniMarinacciRustichini:2006b` define the domain of their cost function to be probabilities over the state space. In our analysis, the state space is $\Theta$, which means that their application of variational preferences gives rise to the robust Bayesian approach in [section](sec:notknowingprior1).
````

(misspecified_likelihood)=
## A misspecified likelihood function

We now propose a generalization of the previous approach by starting from a parameterized family of probabilities $ \tau(w \mid \theta)$ and prior probability measure $\pi$. Typically, a family of parameterized family of probability models is specified so that each model is absolutely continuous with respect to an underlying measure, a condition required to apply likelihood-based methods. Consider relative densities ${\hat m}$ that for each $\theta$ have been rescaled so that

```{math}
:label: rel_likelihood_set
\int_W {\hat m}(w \mid \theta) \tau (w \mid \theta) d\upsilon(w) = 1.
```

To acknowledge misspecification of a model implied by parameter $\theta$, let  ${\hat m}(w \mid  \theta)$  represent an "unstructured" relative perturbation with a parameterized family of densities: 

```{math}
{\hat \ell}( w \mid \theta) = {\hat m}(w \mid  \theta) \tau( w \mid \theta)
```

where ${\hat \ell}(\cdot \mid \theta) \in {\mathcal L}$ for each $\theta \in \Theta.$ With this in mind, let ${\widehat {\mathcal M}}$ be the space of admissible relative densities ${\hat m}(w \mid \theta)$ associated with model $\theta$ for each $\theta \in \Theta$. The pair $({\hat m}, \theta)$ implies a probability distribution represented as

```{math}
{\hat m}(w \mid \theta) \tau(w \mid \theta) d\upsilon(w)
```

over $W$ conditioned on $\theta$. When ${\hat m}$ is not identically one, we view this as a misspecified likelihood function. Uncertainty about the nature of this misspecification induces corresponding uncertainty in the induced distribution,  or the lottery in the language of decision theory.

Preferences that acknowledge this form of model misspecification are ordered by solutions to

```{math}
:label: eqn:minproblem108
\min_{{\hat m}  \in {\widehat {\mathcal M}}} \left( \int_W u[\gamma(w)] {\hat m}(w \mid \theta) \tau(w \mid \theta) d\upsilon(w)   + \xi \int_W \phi[{\hat m}(w \mid \theta)] \tau(w \mid \theta) d \upsilon(w) \right) d\pi_o(\theta)   ,
```

````{prf:remark}
Please remember that we have left dependencies of $\tau$ and $\gamma$ on $\delta$ implicit. Consequently, constraint {eq}`rel_likelihood_set` holds for each $\delta\in\Delta$, where $\tau$ depends implicitly on $\delta.$
````

````{prf:remark}
Another approach would be to use the baseline prior to construct:
```{math}
\int_\Theta \tau(w \mid \theta) d \pi_o(\theta)
```
and treat this as the baseline model of [section](sec:robustnesstypesof); this would correspond to a predictive distribution provided that learning is not formally incorporated into the analysis or that the "prior" $d \pi_o$ has already conditioned on what has been learned from available data.[^Chamberlain2020]

[^Chamberlain2020]: See {cite}`Chamberlain:2020` for a discussion of robustness relative to a predictive distribution.
````

(sec:robustnesstypesof)=
## Robustness reconsidered

It is useful to compare two approaches to robustness that we have taken. The subsection [Robust preferences](sec:robex1) decision maker starts with a baseline prior over parameter vectors and considers consequences of misspecifying that prior. This decision maker takes as given the parameterized family of densities $\tau( w \mid \theta)$ for $\theta \in \Theta$. In contrast, the [section](sec:notknowingmodel) decision maker searches over the entire space $\widehat {\mathcal M}$, subject to a penalty on a statistical divergence from a baseline parameterized family of models. This decision maker considers only the baseline prior distribution.

Our setup allows the parameter space to be infinite dimensional. Consider a prior $\pi_o$ that is consistent with a Bayesian approach to "nonparametric" estimation and inference. Since $\tau(\cdot \mid \theta)$ can be viewed as a mapping from $\Theta$ into ${\mathcal L}$, a prior distribution $\pi_o$ over $\Theta$ implies a corresponding distribution over ${\mathcal L}$. This procedure necessarily assigns prior probability zero to a substantial portion of the space ${\mathcal L}$. Specifying a prior over the infinite dimensional space ${\mathcal L}$ brings challenges associated with all nonparametric methods, including "nonparametric Bayesian" methods that must assign probability one to what is called a "meager set."[^meagerSet] A meager set is defined topologically as a countable union of nowhere dense sets and is arguably small within an infinite-dimensional space. This conclusion carries over to situations with families of priors that are absolutely continuous with respect to a baseline prior, as we have here. To us, prior robustness of this form is interesting, although it is distinct from robustness to potential likelihood misspecifications. Indeed, the [section](sec:notknowingmodel) decision maker who is concerned about model misspecification does not restrict himself to priors that are absolutely continuous with respect to a baseline prior because doing so would exclude many probability distributions he is concerned about.

[^meagerSet]: \citet{sims2010understanding} critically surveys an extensive statistical literature on this issue. Foundational papers are \citet{Freedman_1963}, \citet{SimsAnnals71}, and \citet{Diaconis_Freedman_1986}.

The distinct ways in which the [section](sec:nonknowingprior) and [section](sec:notknowingmodel) formulations use statistical discrepancies lead to substantial differences in the associated variational preferences, namely, representation {eq}`eqn:minproblem106` or {eq}`eqn:minproblem107` for the [section](sec:nonknowingprior) way of prior ambiguity and representation {eq}`eqn:minproblem108` way of ambiguity about the parameterized family of densities, $\tau( \cdot \mid \theta)$.

(sec:example)=
## Two examples

It is instructive to apply the distinct approaches of [sections](sec:nonknowingprior) and [section](sec:notknowingmodel) to simple examples. The first example gives a simple illustration of preference inputs into robust control problems, and the second one explores a forecasting problem.

(sec:robex1)=
### Robust preferences

Assume the following constituents:
- Baseline model is $\tau_o(w) \sim \text{Normal} (\mu_o, \sigma_o^2)$
- Alternative structured models $\tau(w \mid \theta_i) \sim \text{Normal} (\mu_i, \sigma_i^2), i = 1, \ldots, k$, where potential parameter values (states) are $\theta_i = (\mu_i, \sigma_i)$ and parameter space $\Theta = \{\theta_i : i=1,2, \ldots, k\}.$. The baseline model can be one of these $k$ models.
- Baseline prior over structured models is a uniform distribution $\pi_o(\theta_i) = \frac{1}{k}, i = 1, \ldots, k.$
- Prize is the induced distribution of $c(w) = \gamma(w).$
- Utility function is $u[c(w)] = \log [c(w)]$, where $c(w)$ is consumption
- Prize rule is $\gamma(w) = \exp(\gamma_0 + \gamma_1 w)$

To obtain an alternative prior $\pi_i$ for $i=1, \ldots, k$, we set $n_i = k \pi_i$ so that the product of $n_i$ times the baseline prior is:
```{math}
\frac{n_i}{k} = \pi_i.
```
The expected utility conditioned on parameter vector $\theta_i$ is
```{math}
\int_{W} u[\exp (\gamma_0 + \gamma_1 w)] \tau(w | \theta)d\upsilon(w) = \gamma_0 + \gamma_1 \mu_i,
```
and a statistical divergence applied to alternative priors is
```{math}
{\frac{1}{k}} \sum_{i=1}^k \phi(k \pi_i) .
```
A [subsection](sec:notknowingprior1) decision maker with variational preferences orders prize rules $\gamma(w) = \exp(\gamma_0 + \gamma_1 w)$ according to
```{math}
\min_{\pi_i\ge 0, \sum_{i=1}^k \pi_i = 1} \gamma_0 + \gamma_1 \sum_{i=1}^k \pi_i \mu_i +
{\frac{\xi}{k}}  \sum_{i=1}^k \phi(k \pi_i) .
```
For a relative entropy divergence, prize rules are ordered by
```{math}
-\xi \log \sum_{i=1}^k \left({\frac{1}{k}}\right)  \exp\left[-{\frac{1}{\xi}} \left(\gamma_0 + \gamma_1 \mu_i \right) \right]
= \gamma_0 -\xi \log \sum_{i=1}^k \left({\frac{1}{k}}\right) \exp\left(-{\frac{\gamma_1\mu_i}{\xi}}\right),
```
and the associated minimizing $\pi_i$ is
```{math}
{\frac{\exp\left(-{\frac{\gamma_1 \mu_i}{\xi}}\right)}{\sum_{i=1}^k \exp\left(-{\frac{\gamma_1 \mu_i}{\xi}}\right)}}.
```
A [subsection](sec:nonknowingprior2) decision maker, in effect, chooses the multiplier $\xi$ to hit a relative entropy constraint on the prior.

A criterion that expresses robustness to prior misspecification with a relative entropy divergence ranks prizes as either
```{math}
-\xi \log \sum_{i=1}^k \left(\frac{1}{k}\right) \exp\left[-\frac{1}{\xi} \left(\gamma_0 + \gamma_1 \mu_i \right)\right],
```
or,
```{math}
\max _\xi -\xi \log \sum_{i=1}^k \left({\frac{1}{k}}\right) \exp\left[-\frac{1}{\xi} \left(\gamma_0 + \gamma_1 \mu_i \right)\right] - \xi \kappa
```

When we use relative entropy as a statistical divergence, variational preferences for a [subsection](sec:notknowingmodel) decision maker are ordered by
```{math}
\gamma_0 + \gamma_1 \mu_0 - {\frac{1}{2\xi}} (\sigma_0 \gamma_1)^2
```
Larger values of the positive scalar $\xi$ call for smaller adjustments $- {\frac{1}{2\xi}} (\sigma_o \gamma_1)^2$ of expected utility $\gamma_0 + \gamma_1 \mu_o$ for concerns about misspecification of $\tau_o(w)d\upsilon(w)$.

### Robust forecasting

Consistent with Example {prf:ref}`remark:Ferguson`, partition
```{math}
w = \begin{bmatrix} w_1 \\ w_2 \end{bmatrix},
```
where $w_1$ is scalar outcome of a variable to be forecast and $w_2$ constitutes data underlying a forecast. Assume that:
- The baseline model is $\tau_o(w)$. 
- Alternative structured models are $\tau(w \mid \theta)$ for a parameter space $\Theta = \{\theta_i : i=1,2, \ldots, k\}.$ The baseline model can be one of the $k$ models.
- The baseline prior over structured models is a uniform distribution $\pi_o(\theta_i) = \frac{1}{k}, i = 1, \ldots, k.$
- The prize is the induced distribution of the forecast error $w_1 - \delta(w_2)$, where $\delta$ is the forecast rule.
- The utility function is $-[w_1 - \delta(w_2)]^2.$
- The prize rule is $\gamma_\delta(w) = w_1 - \delta(w_2).$

(sec:hybrid)=
## Hybrid models

We now use components described above as inputs into a representation of preferences that includes uncertainty about a prior to put over structured models as well as concerns about possible misspecifications of those structured models. We use probability perturbations in the form of alternative relative densities in ${\widehat {\mathcal M}}$ to capture uncertainty about models and probability perturbations in the form of alternative relative densities ${\mathcal N}$ to capture uncertainty about a prior over models.

Let $\pi_o(\theta)$ is a baseline prior over $\theta$. To conduct a prior robustness analysis, consider alternative priors

```{math}
d \pi(\theta) = n(\theta) d \pi_o(\theta)
```

for $n \in {\mathcal N}.$

Consider relative densities ${\hat m}$ that for each $\theta$ have been rescaled so that

```{math}
\int_W {\hat m}(w \mid \theta) \tau(w \mid \theta) d\upsilon(w) = 1.
```

To acknowledge misspecification of a model implied by parameter $\theta$, let ${\hat m}(w \mid \theta)$ represent an "unstructured" perturbation of that model. With this in mind, let ${\widehat {\mathcal M}}$ be the space of admissible relative densities ${\hat m}(w \mid \theta)$ associated with model $\theta$ for each $\theta \in \Theta$. We then consider a composite parameter $({\hat m}, \theta)$ for ${\hat m} \in {\widehat {\mathcal M}}$ and $\theta \in \Theta.$ The composite parameter $({\hat m}, \theta)$ implies a distribution ${\hat m}(w \mid \theta) \ell(w \mid \theta) \tau(w \mid \theta)d\upsilon(w)$ over $W$ conditioned on $\theta$.

To measure a statistical discrepancy that comes from applying $\hat m$ to the density $\ell$ of $w$ conditioned on $\theta$ and by applying $n$ to the baseline prior over $\theta$, we first acknowledge possible misspecification of each of the $\theta$ models by computing:

```{math}
{\mathbb T}_1[\gamma]( \theta)  = \min_{{\hat m} \in {\widehat {\mathcal M}}} \int_W \left(  u[\gamma(w)]{\hat m}(w  \mid \theta)  + \xi_1 \phi_1\left[ {\hat m}(w \mid \theta)  \right] \right) \ell( w \mid \theta)  \tau(w \mid \theta)d\upsilon(w)
```

The ${\mathbb T}_1$ operator maps prize rules $\gamma$ into functions of $\theta$. We use this for both hybrid approaches.

(hybrid_one)=
### First hybrid model

We can rank alternative prize rules $\gamma$ by including the following adjustment for possible misspecification of the baseline prior $\pi_o$:

```{math}
{\mathbb T}_2 \circ {\mathbb T}_1[\gamma] = \min_{ n \in {\mathcal N}} \int_\Theta \left({\mathbb T}_1[ \gamma]( \theta) n(\theta) + \xi_2 \phi_2[n(\theta)] \right) d \pi_o(\theta).
```

Here $\phi_1$ and $\phi_2$ are possibly distinct convex functions with properties like the ones that we imposed on $\phi$ in [section](sec:statdivergasc).

Such a two-step adjustment for possible misspecification leads to an implied one-step variational representation with a composite divergence that we can define in the following way. For ${\hat m} \in \widehat{\mathcal M}$ and $n \in {\mathcal N}$, form a composite scaled statistical discrepancy

```{math}
:label: double_penalty
{\widehat  D}( {\hat m}, n \mid \tau, \pi_o) = \xi_1 \int_\Theta \left(  \int_W \phi_1\left[ {\hat m}(w \mid \theta) \right] \tau(w \mid \theta)  d\upsilon(\theta) \right) n(\theta) d\pi_o(\theta)  + \xi_2 \int_\Theta \phi_2 \left[ n(\theta) \right] d \pi_o(\theta)
```

for $\xi_1 > 0, \xi_2 > 0$. Then variational preferences are ordered by

```{math}
\min_{{\hat m} \in {\widehat {\mathcal M}}, n \in {\mathcal N} } \int_\Theta \left(\int_W u[\gamma(w)] {\hat m}(w \mid \theta)  \tau( w \mid \theta) d\upsilon(w) \right)  n(\theta) d\pi_o(\theta) + {\widehat D}({\hat m}, n \mid \tau, \pi_o)
```

In [Appendix](app:convexity) we establish that divergence {eq}`double_penalty` is convex over the family of probability measures that concerns the decision maker.

````{prf:remark}
As noted earlier, {cite}`Cerreiaetal:2013` posit a state space that includes parameters but also can include what we call repercussions. Thus, think of the state as the pair $(w,\theta)$. In this setting, one could apply a statistical divergence to a joint distribution over possible realizations of $(w,\theta)$. Since the joint distribution can be factored into the product of a distribution over $W$ conditioned on $\theta$ and a marginal distribution over $\Theta$, such an approach can capture robustness in the specification of both $\tau$ and $\pi_o$, albeit in a very specific way. For instance, for the relative entropy divergence, this results in the joint divergence measure:

```{math}
{\widehat D}( {\hat m}, n \mid \tau, \pi_o ) =  \xi_1 \int_\Theta \left[  \int_W {\hat m}(w \mid \theta) \log{\hat m}(w \mid \theta)   \tau(w \mid \theta) d\upsilon(w) \right] n(\theta) d\pi_o(\theta)  + \xi_2
\int_\Theta n(\theta) \log n(\theta)d \pi_o(\theta)
```
for $\xi_1 = \xi_2$.

In earlier work, we have demonstrated important limits to such an approach in dynamic settings.[^limits-dynamic-settings] As we have shown here, we find both robustness to model misspecification and robustness to prior specification to be interesting in their own rights and see little reason to group them into a single $\phi$ divergence.
````

[^limits-dynamic-settings]: See {cite}`HansenSargent:2007,HansenSargent:2011` and {cite}`HansenMiao:2018`.

(hybrid_two)=
## Second hybrid model
As an alternative to the [approach](hybrid_one) approach, we could instead constrain the set of priors to satisfy:

```{math}
:label: kappa_constraint
\int_\Theta \phi_2[n(\theta)]  d\pi_o(\theta) \le \kappa
```
so that a decision maker's preferences over prize rules $\gamma$ would be ordered by:

```{math}
:label: CVMM103
\min_{n \in {\mathcal N}} \int_\Theta {\mathbb T}_1\left[ \gamma\right](\theta) n(\theta) d\pi_o(\theta),
```
where minimization is subject to {eq}`kappa_constraint`.

As in {cite}`Cerreiaetal:2021`, preferences ordered by {eq}`CVMM103` subject to constraint {eq}`kappa_constraint` can be thought of as using a divergence between a potentially misspecified probability distribution and a set of predictive distributions that have been constructed from priors over a parameterized family of probability densities within the constrained set $\Theta`.[^cerreia-2021-axiomatic]

Notice how the first term in discrepancy measure {eq}`double_penalty` uses a prior $n d\pi_o$ to construct a weighted averaged over $\theta \in \Theta$ of the following conditioned-on-$\theta$ misspecification measure
$$
\xi_1 \left(  \int_W \phi_1\left[ {\hat m}(w \mid \theta) \right] \tau(w \mid \theta) d \upsilon(w) \right).
$$
The objective in problem {eq}`CVMM103` is to make the divergence between a given distribution and each of the parameterized probability models small on average by minimizing over how to weight divergence measures indexed by $\theta$ subject to the constraint that $\pi \in \Pi$.[^structured-models] Equivalently, in place of {eq}`double_penalty`, this approach uses cost function
$$
{\widetilde D} ( {\hat m} \mid \tau, \pi_o) = \xi_1 \min_{n \in {\mathcal N}} \int_{\Theta} \left(  \int_W \phi_1\left[ {\hat m}(w \mid \theta) \right] d\ell(w \mid \theta) \right) n(\theta) d\pi_o(\theta) .
$$

[^cerreia-2021-axiomatic]: {cite}`Cerreiaetal:2021` provide an axiomatic justification of set-based divergences as a way to capture model misspecification within a {cite}`Gilboaetal:2010` setup with multiple models.
[^structured-models]: By emphasizing a family of structured models, this set-divergence concept differs from an alternative that could be constructed in terms of an implied family of predictive distributions.

````{prf:remark}
It is possible to simplify computations by using dual versions of the hybrid approaches delineated in subsections [hybrid approaches one](hybrid_one) and [hybrid approaches two](hybrid_two). Such formulations closely parallel those described in our discussions of robust prior analysis and potential model misspecification in remarks {prf:ref}`rem:dual1`, {prf:ref}`rem:dual2`, and {prf:ref}`rem:dual3`.
````

(sec:dynext)=
## Dynamic extension

Although a complete treatment of dynamics deserves its own paper, here we describe briefly how to extend the familiar recursive utility specification of {cite}`KrepsPorteus:1978` and {cite}`EpsteinZin:1989` to accommodate our two robustness concerns to an intertemporal environment. We accomplish this by using conditional counterparts to the preceding analysis to explore consequences of mis-specifying Markov transition dynamics and prior distributions over unknown parameters. The resulting preferences have a recursive structure. There is an inherent tension between dynamic consistency and statistical consistency in these preferences that we discuss elsewhere ({cite}`HansenSargent:2022decision`).  

### A deterministic warm up

We represent preferences using recursions that apply to continuation values. Abstracting from uncertainty, a commonly used intertemporal preference specification is captured by the value recursion:
```{math}
V_t = \left[(1 - \beta) \left(C_t\right)^{1-\rho} + \beta \left( V_{t+1} \right)^{1 - \rho} \right]^{\frac 1 {1-\rho}}
```
for $0<\beta<1$ and $\rho> 0$. $V_t$ is the date $t$ continuation value and $C_t$ is date $t$ consumption. The parameter $\beta$ governs discounting and the parameter $\rho$ is the reciprocal of the intertemporal elasticity of substitution. Applying the recursion over an infinite horizon leads to the following expression for the continuation value:
```{math}
V_t = \left[ (1-\beta) \sum_{j=0}^\infty \beta^j \left(C_{t+j} \right)^{1-\rho}\right]^{\frac 1 {1-\rho}}
```
Since the logarithmic transformation is increasing, we can use the following recursion in the logarithm ${\widehat V}_t$ of the continuation value to represent preferences:
```{math}
{\widehat V}_t = \frac 1 {1-\rho} \log 
\left[(1 - \beta) \exp\left[(1-\rho){\widehat C}_t\right] + \beta \exp \left[(1-\rho){\widehat V}_{t+1} \right] \right] 
```
where ${\widehat C}_t$ is the logarithm of consumption.  

### Introducing uncertainty

Let ${\mathfrak A}_t$ denote a sigma algebra capturing information available to the decision maker at date $t$. Think of the repercussion $W_{t+1}$ as generating new information relative to ${\mathfrak A}_t$ that is pertinent for constructing ${\mathfrak A}_{t+1}$. Think of the continuation value, ${\widehat V}_{t+1}$ as the counterpart to a prize that can depend on a repercussion vector $W_{t+1}.$ A continuation value ${\widehat V}_{t+1}$ is constrained to be measurable with respect to ${\mathfrak A}_{t+1}$. We explore model misspecification by using nonnegative random variables $M_{t+1}$ that are ${\mathfrak A}_{t+1}$ measurable and satisfy ${\mathbb E} \left( M_{t+1} \mid {\mathfrak A}_t, \theta\right) = 1.$ We explore prior/posterior misspecification using nonnegative random variables $N_{t}$ that are measurable with respect ${\mathfrak A}_t$ augmented by knowledge of $\theta$ and satisfy ${\mathbb E}\left( N_t \mid {\mathfrak A}_t \right) = 1.$

To accommodate robustness concerns in decision making, define preferences with three recursions for updating the continuation value  
```{math}
:label: three_recur
\begin{align}
{\widehat V}_t &= \frac 1 {1-\rho} \log 
\left[(1 - \beta) \exp\left[(1-\rho){\widehat C}_t\right] + \beta \exp \left[(1-\rho){\overline R}_t \right] \right] \\
{\widehat R}_t & = \min_{M_{t+1}\ge 0, {\mathbb E} \left(M_{t+1} \vert {\mathfrak A}_t, \theta \right) = 1}  
{\mathbb E} \left[ {M_{t+1} \widehat V}_{t+1} + \xi_1 \phi_m \left(M_{t+1} \right) \mid {\mathfrak A}_t, \theta \right]  \\
{\overline R}_t & = \min_{N_{t}\ge 0, {\mathbb E} \left(N_{t} \vert {\mathfrak A}_t \right) = 1}  
{\mathbb E} \left[ N_{t} {\widehat R}_{t} + \xi_2 \phi_n \left(N_{t} \right) \mid {\mathfrak A}_t \right] 
\end{align}
```
The second and third recursions provide a dynamic counterpart to the approach in [section hybrid_one](sec:hybrid_one). Replacing the third recursion in {eq}`three_recur` with a constrained counterpart gives a dynamic counterpart to the approach in [section hybrid_two](sec:hybrid_two).[^footnoteAltApproach]

[^footnoteAltApproach]: See {cite}`HansenSargent:2020, HansenSargent:2022decision` for elaboration and application of this alternative approach.
### Shadow valuation

Following {cite}`HansenRichard:1987` and others, we can use stochastic discount factors to value assets having an uncertain one-period ahead payoff. We deduce shadow values by computing a one-period intertemporal marginal rate of substitution. Of particular interest to us are contributions that our model-misspecification operator
${\widehat R}_t$ and our prior-robustness operator ${\overline R}_t` make to this shadow value.   

A contribution to the shadow value that comes from the first recursion in {eq}`three_recur` looks at marginal contributions in adjacent time periods.  
Date $t$ marginal contributions of $C_t$ and ${\overline R}_t$ to the current period continuation value are:
```{math}
\begin{align*}
MC_t & = 
(1-\beta) \exp\left[ (\rho - 1) {\widehat V}_t \right] 
\left(C_t\right)^{-\rho} \\
M{\overline R}_t & = \beta \exp\left[ (\rho - 1) {\widehat V}_t \right] \exp\left[ (1-\rho) {\overline R}_t \right] .
\end{align*}
```
Since our aim is to infer the one-period intertemporal marginal rate of substitution, we look across adjacent time periods using consumption at each date as a numeraire:
```{math}
\frac {{MC_{t+1} M{\overline R}_t }}{MC_t} = \beta \left( \frac {C_{t+1}}{C_t} \right)^{-\rho} \exp\left[(\rho - 1) \left( {\widehat V}_{t+1} - {\overline R}_t \right) \right]. 
```
This would give the deterministic intertemporal marginal rate of substitution if we were to substitute ${\widehat V}_{t+1}$ for ${\overline R}_{t}$ in this expression.     

For the uncertainty adjustments, we deduce the marginal contributions by applying the Envelope Theorem to the minimization problems in the second and third recursions in {eq}`three_recur`: 
- $M{\widehat V}_{t+1} = M_{t+1}^*$
- $M {\widehat R}_t = N_{t}^*$

Thus, the minimizing changes in probabilities contribute directly to the shadow valuation. The resulting increment to a stochastic discount factor process is:
```{math}
\frac {S_{t+1}}{S_t} = \beta \left( \frac {C_{t+1}}{C_t} \right)^{-\rho} \exp\left[(\rho - 1) \left( {\widehat V}_{t+1} - {\overline R}_{t} \right) \right]M_{t+1}^*N_{t}^*
```
where 
- $M_{t+1}^*$ adjusts for possible model misspecification
- $N_{t}^*$ adjusts for possible prior misspecification 

(sec:uncert_quant)=
## An approach to uncertainty quantification

[Subsection](sec:hybrid) posed a minimum problem that comes from variational preferences with a two-parameter cost function that we constructed from two statistical divergences. Along with a robust prize rule, the minimum problem produces a worst-case probability distribution that rationalizes that prize rule. Strictly speaking, the decision theory tells us that particular values of cost function parameters $(\xi_1, \xi_2)$ express a decision maker's concerns about uncertainty, broadly conceived. In the spirit of {cite}`Good:1952`, it can be enlightening to study how worst-case distributions depend on $(\xi_1, \xi_2)$. The concluding paragraph of {cite}`Chamberlain:2020` recommends exploring sensitivities with respect to a likelihood and with respect to a prior. Sensitivity of worst-case distributions to $(\xi_1, \xi_2)$ provides evidence about the forms of subjective uncertainty and potential model misspecification that *should* be of most concern. That can provide decision makers and outside analysts better understandings of the consequences of uncertainty aversion.

Motivated partly by a robust Bayesian approach, we have used decision theory to suggest a new approach to uncertainty quantification. By varying the aversion parameters $(\xi_1, \xi_2),$ we can trace out two-dimensional representations of prize rules and worst-case probabilities. A representation of worst-case probabilities includes both worst-case priors and a worst-case alteration to each member of a parametric family of models. A decision maker can explore alternative choices and associated expected utilities by studying how $(\xi_1, \xi_2)$ trace out a two-dimensional set of worst-case probabilities. In this way, we reduce potentially high-dimensional subjective uncertainties to a two-dimensional collection of alternative probability specifications that should most concern a decision maker along with accompanying robust prize rules for responding to those uncertainties. 

(sec:statlearning)=
## Relation to statistical learning

We briefly compare our approach to related analyses coming from statistical learning theory and, in particular, PAC (probably approximately correct) Bayesian analysis. See {cite}`Guedj:2019` for a recent survey of PAC Bayesian methods and see {cite}`McAllester:1999` and {cite}`Cantani:2007`, among others, for fundamental contributions. While their formulations of a decision problem differ from ours, there are intriguing connections.

To understand some of the connections, partition a random vector, $Y$, with realization $y$ as 
```{math}
Y' :=  \begin{bmatrix} {Y_1}' & {Y_2}' &  {Y_2}' & ... & {Y_K}' \end{bmatrix},
```
and regard it as  "training data" for a machine learning method.   
For an objective function, construct an "empirical risk" criterion:
```{math}
{\widehat \Phi }(Y, \theta)  := {\frac 1 K} \sum_{k=1}^K  \Phi (Y_k, \theta) .
```
The object of interest $\theta$ can be an element of a collection of functions.  
Define the population counterpart to ${\widehat \Phi }(Y, \theta)$ as  
```{math}
{\overline \Phi} (\theta ),
```
the Law of Large Numbers limit of ${\widehat \Phi }(W, \theta)$ as $K \rightarrow \infty$. Suppose that an idealized target decision solves:
```{math}
\theta^* = \arg \min_{\theta \in {\Theta}}  {\overline{\Phi}} (\theta ).
```
This estimator appears in an extensive literature on $M$ estimation; the idealized optimized decision $\theta^*$ defines a parameter or decision of interest. In this setting, we use decisions and parameters interchangeably, in contrast to our formulation.

A common approach to M estimation is to solve the finite sample analog problem
```{math}
{\hat \theta} = \arg \min_{\theta \in \Theta} {\widehat \Phi }(Y, \theta) .
```
This approach struggles when the space $\Theta$ is “large”, the typical case with machine learning methods that fit flexible functional forms with many parameters. More generally, statistical learning often seeks meaningful worst-case bounds of finite sample approximates to a solution of the population problem.

Motivated by concerns for applications when the space $\Theta$ in standard M estimation is expansive, the PAC Bayesian approach proceeds differently. The approach seeks a probability distribution $\pi$ over the space $\Theta$ given the data $W$, rather than a single value, $\theta$. By analogy to Bayesian methods, such a distribution is referred to as a `generalized posterior.' The approach imposes a baseline prior distribution $\pi_o$ over the space $\Theta$, and considers generalized posteriors in a family:
```{math}
d \pi(\theta) = n(\theta) d\pi_o(\theta)
```
for $\int n d\pi_o = 1.$

Instead of solving the finite sample M estimation problem, consider a family of problems:

```{math}
:label: PAC-problem
min_{n, \int n d\pi_o = 1}   \int_{{\Theta}} {\widehat \Phi }(Y, \theta) d \pi(\theta) + \xi \int_{\Theta} \log n(\theta) n(\theta) d \pi_o(\theta) 
```
indexed by $\xi$. Applying the same mathematics we have used in previous sections, minimization brings exponential tilting:
```{math}
n^*(\theta) = \frac{ \exp \left[ - \frac 1 \xi  {\widehat \Phi }(Y, \theta)\right]}
{\int_\Theta  \exp \left[ - \frac 1 \xi  {\widehat \Phi }(Y, {\tilde \theta})\right] d \pi_o({\tilde \theta}) }
```
The Bayesian PAC uses this minimizer to construct an approximation to a minimizer of the underlying (infeasible) population problem. 

Problem {eq}`PAC-problem` provides a way to incorporate probabilistic restrictions into the M estimation problem. There are some interesting special cases. When $ {\widehat \Phi }(W, \cdot)$ is the negative of the log likelihood function and $\xi = 1$, we are led to a standard calculation of a Bayesian posterior. {cite}`zhang:2006`, {cite}`grunwald:2011` and others propose and defend a `safe Bayesian framework' by exploring other values of $\xi > 1$ based on robustness considerations. The Bayesian PAC approach studies alternative specifications of ${\widehat \Phi }(W, \cdot)$ based on more general loss functions. As in M-estimation more generally, this construction may embed some robustness concerns. When $\xi$ tends to infinity the generalized posterior collapses to the prior distribution. When $\xi$ tends to zero the prior becomes inconsequential, and the generalized posterior collapses to the solution the finite-sample M estimation solution. The penalty parameter $\xi$ governs a tradeoff between the importance of the objective ${\widehat \Phi }(Y, \cdot)$ and the baseline prior $\pi_o$. The PAC-Bayesian literature discusses extensively the role of the parameter in approximation.

While our approach shares much mathematical structure with PAC-Bayesian methods, it differs in ways that are significant for applications. The M estimation formulation ties its decision problem directly to an underlying unknown "parameter" and contains no counterpart to the maximization steps that we use to represent uncertainty aversions. Furthermore, the PAC-Bayesian Problem {eq}`PAC-problem` conditions on $W$ and focuses exclusively on uncertainties about unknown states or parameters. Also, PAC-Bayesian methods use Problem {eq}`PAC-problem` as a device to approximate a solution to an infeasible population problem, which is not a component of our analysis.

To elaborate more on the differences between PAC-Bayesian methods and our approach, we study decision problems in which parameters are not the objects of ultimate interest but instead are just intermediate "means to ends" of constructing decision rules for making choices of economics quantities that are robust to misspecifications. Rather than replacing a log-likelihood function with an M estimation objective and possibly down-weighting its importance, we introduce potential likelihood misspecifications explicitly; we also formally acknowledge possible misspecification of priors. By appropriately adjusting the divergence “cost” structure, our approach allows us to explore tradeoffs between concerns about misspecifications of likelihoods, on the one hand, and priors, on the other hand. 

(sec:conclude)=
## Concluding remarks

Except for our brief section [excursion](sec:dynext), we have confined ourselves to a "static" setting that allowed us to apply and extend a framework created by {cite}`MaccheroniMarinacciRustichini:2006b` to distinguish ambiguity about a prior and from concerns for misspecifications of likelihood functions. In doing this we reinterpret objects that appear in the {cite}`Anscombe_Aumann` formulation to represent both types of doubts, distinct specification concerns that are familiar to applied statisticians. We intend the present paper as a prolegomenon to a sequel in which we shall extend and reinterpret the dynamic variational preferences of {cite}`MaccheroniMarinacciRustichini:2006`. That dynamic formulation will connect to a dynamic measure of statistical divergence based on relative entropy and the recursive preferences of {cite}`KrepsPorteus:1978` and {cite}`EpsteinZin:1989`. While the issues studied here will arise in that framework, additional ones such as dynamic consistency and appropriate choices of state variables for recursive formulations of preferences also appear.[^dynamicAssetPricingNote]

[^dynamicAssetPricingNote]: To apply quantum methods to dynamic asset pricing models, {cite}`GhyselsMorgan:2023` deploy extensions of the formulation developed here. An interesting notion of a "state" again comes into play.

# Appendix

(app:convesity)=
## Convexity of composite divergence 

To verify convexity of {eq}`double_penalty`, consider two joint probability measures on $W \times \Theta$:
```{math}
&{\hat m}_0(w \mid \theta) \tau(w \mid \theta) d \upsilon(w) n_0(\theta) d \pi_o(\theta) \\
&{\hat m}_1(w \mid \theta) \tau(w \mid \theta) d \upsilon(w) n_1(\theta) d\pi_o(\theta) .
```
A convex combination of these two probability measures is itself a probability measure. Use weights $1 - \alpha$ and $\alpha$ to construct a convex combination and then factor it in the following way. First, compute the marginal probability distribution for $\theta$ expressed as $n_\alpha(\theta) d\pi_o(\theta)$:
```{math}
n_\alpha (\theta) = (1- \alpha) n_0(\theta) +\alpha n_1(\theta) .
```
By the convexity of $\phi_2$, it follows that
```{math}
:label: first_inequality
\phi_2[n_\alpha (\theta) ] \le (1 - \alpha) \phi_2 [n_0(\theta)] + \alpha \phi_2[n_1(\theta)] .
```
Next note that
```{math}
{\hat m}_\alpha(w \mid \theta) = & \left[{\frac {(1-\alpha) n_0(\theta)}  {(1-\alpha) n_0(\theta) + \alpha n_1(\theta)}} \right]
{\hat m}_0(w \mid \theta) \\ & + \left[{\frac { \alpha n_1(\theta)}  {(1-\alpha) n_0(\theta) + \alpha n_1(\theta)}} \right]
{\hat m}_1(w \mid \theta) .
```
By the convexity of $\phi_1$
```{math}
\phi_1[{\hat m}_\alpha(w \mid \theta)] \le & \left[{\frac {(1-\alpha) n_0(\theta)}  {(1-\alpha) n_0(\theta) + \alpha n_1(\theta)}} \right] \phi_1 [{\hat m}_0 (w \mid \theta) ] \\ & + \left[{\frac { \alpha n_1(\theta)}  {(1-\alpha) n_0(\theta) + \alpha n_1(\theta)}} \right]
\phi_1[  {\hat m}_1(w \mid \theta) ].
```
Thus,
```{math}
:label: second_inequality
\phi_1[{\hat m}_\alpha(w \mid \theta)]n_\alpha(\theta) \le (1-\alpha) n_0(\theta) \phi_1 [{\hat m}_0 (w \mid \theta) ] +  \alpha n_1(\theta)\phi_1[  {\hat m}_1(w \mid \theta) ] .
```
Multiply {eq}`second_inequality` by $\xi_1$ and {eq}`first_inequality` by $\xi_2$, add the resulting two terms, and integrate with respect to $\tau(w \mid \theta) d \upsilon(w)  d \pi_o(\theta)$ to verify that divergence {eq}`double_penalty` is indeed convex in probability measures that concern the decision maker.

