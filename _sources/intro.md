# Quant Macro Finance
Quant Macro Finance (QMF) is a series of self-contained online texts intended to accompany Professor Lars Peter Hansen's lectures on Asset Pricing I.

```{image} images/Cartoon.png
:alt: fishy
:class: bg-primary mb-1
:align: center
```

<!-- QMF provides students with the pedagogical and coding resources required to solve dynamic stochastic general equilibrium models under uncertainty and to compute corresponding:
<br>
- <i>Stationary densities</i>, the distribution of variables of interest under stationarity, </li>
- <i>Shock-exposure elasticities</i>, which capture the exposure of macroeconomic processes to shocks over different investment horizons, and </li>
- <i>Shock-price elasticities</i>, which capture changes in risk compensations associated with shock exposures over different investment horizons</li>
<br> -->
First, we explain the theoretical framework for solving nonlinear DSGE models with uncertainty and computing shock elasticites.

1. [Recursive Utility](theory/Recursive.md)
2. [Uncertain Expansion](theory/uncertainexpansion.ipynb)
3. [Shock Elasticities](theory/shockelasticity.ipynb)

Second, we introduce toolboxes used to approximate solutions and compute shock elasticities in discrete time.

4. [Quick Guide to Expansion](discrete_local_approximation/quickguide.ipynb)
5. [Computing Elasticities in Discrete Time]

Then, we provide extensions via continuous time, two-capital models and financial intermediaries. 

6. [Two Capital Model]
7. [Computing Elasticities in Continuous TIme]

<br>
<p>The series is divided into <i>discrete and local</i> and <i>continuous and global</i> methods. Local methods are based on perturbing around the steady state in order to approximate policy and valuation functions in a neighbourhood around the steady state, whereas global methods allow us to compute these functions across the domain.</p>
<br><br>
You can find Professor Hansen's lecture notes <a href = "https://larspeterhansen.org/class-notes/">here</a>. 