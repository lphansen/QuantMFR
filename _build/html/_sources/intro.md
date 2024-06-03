# Quant Macro Finance
Quant Macro Finance (QMF) is a series of online texts intended to accompany *Risk, Uncertainty and Value* by Lars Peter Hansen, Thomas Sargent and Jaroslav Borovička. 

```{figure} images/bernoulli_pissarro.jpg
---
name: Bernoulli-pissarro
---
Over three hundred years ago, Jakob Bernoulli extended the use of probability from games of chance to the study of observable social phenomenon including economic outcomes as illustrated by the Pissarro painting of a marketplace in Rouen    Bernoulli derived a version of the Law of Large Numbers to justify empirical measurements of interest.  See Stephen Stigler’s paper ``Soft Questions, Hard Answers: Jacob Bernoulli’s Probability in Historical Context,’’ published in the International Statistical Review in 2014 for a revealing discussion of Bernoulli’s contribution. 
```

<!-- QMF provides students with the pedagogical and coding resources required to solve dynamic stochastic general equilibrium models under uncertainty and to compute corresponding:
<br>
- <i>Stationary densities</i>, the distribution of variables of interest under stationarity, </li>
- <i>Shock-exposure elasticities</i>, which capture the exposure of macroeconomic processes to shocks over different investment horizons, and </li>
- <i>Shock-price elasticities</i>, which capture changes in risk compensations associated with shock exposures over different investment horizons</li>
<br> -->

Currently, a preliminary HTML version of our first 6 chapters are available: 
1. [Stochastic Processes and Law of Large Numbers](book/example_out_c1_v2.md)
2. [Markov Processes](book/example_out_c2_v2.md)
3. [Stationary Increments](book/example_out_c3_v2.md)
4. [Processes with Markovian Increments](book/example_out_c4_v2.md)
5. [Decision, Ambiguity, and Misspecification](book/decision_book_draft.md)
6. [Hidden Markov Models](book/example_out_c5_v2.md)
7. [Likelihoods](book/example_out_c6_v2.md)

Other notes and material that will eventually be incorporated into chapters follow: 

First, we explain the theoretical background and numerical solution of expansion method to solve nonlinear DSGE models with recursive utility function. 

1. [Uncertain Expansion Theory](theory/uncertainexpansion_update.ipynb)
2. [Numerical Solution: Expansion Suite](theory/quickguide_update.ipynb)

Second, we introduce frameworks to solve price and exposure elasticities in both discrete time and continuous time. 

3. [Shock Elasticity: Discrete Time](continuous_global_solution/shockelasticity.ipynb)
4. [Shock Elasticity: Continuous Time](continuous_global_solution/shockelasticitycontinuous.ipynb) 


<br>
<!-- <p>The series is divided into <i>discrete and local</i> and <i>continuous and global</i> methods. Local methods are based on perturbing around the steady state in order to approximate policy and valuation functions in a neighbourhood around the steady state, whereas global methods allow us to compute these functions across the domain.</p>
<br><br> -->
You can find Professor Hansen's lecture notes and support material for his second year asset pricing class <a href = "https://larspeterhansen.org/class-notes/">here</a>. 