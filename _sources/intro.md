# Quant Macro Finance
Quant Macro Finance (QMF) is a series of online texts intended to accompany *Risk, Uncertainty and Value* by Thomas Sargent and Lars Peter Hansen. 

```{image} images/Cartoon.png
:alt: fishy
:class: bg-primary mb-1
:width: 700px
:align: center
```

<!-- QMF provides students with the pedagogical and coding resources required to solve dynamic stochastic general equilibrium models under uncertainty and to compute corresponding:
<br>
- <i>Stationary densities</i>, the distribution of variables of interest under stationarity, </li>
- <i>Shock-exposure elasticities</i>, which capture the exposure of macroeconomic processes to shocks over different investment horizons, and </li>
- <i>Shock-price elasticities</i>, which capture changes in risk compensations associated with shock exposures over different investment horizons</li>
<br> -->

First, we explain the theoretical background and numerical solution of expansion method to solve nonlinear DSGE models with recursive utility function. 

1. [Uncertain Expansion Theory](theory/uncertainexpansion.ipynb)
2. [Numerical Solution: Expansion Suite](theory/quickguide.ipynb)

Second, we introduce frameworks to solve price and exposure elasticity under discrete time and continuous time settings
4. [Shock Elasticity: Discrete Time](discrete_local_approximation/shockelasticity.ipynb)
5. [Shock Elasticity: Continuous Time](discrete_local_approximation/shockelasticitycontinuous.ipynb) 


<br>
<!-- <p>The series is divided into <i>discrete and local</i> and <i>continuous and global</i> methods. Local methods are based on perturbing around the steady state in order to approximate policy and valuation functions in a neighbourhood around the steady state, whereas global methods allow us to compute these functions across the domain.</p>
<br><br> -->
You can find Professor Hansen's lecture notes <a href = "https://larspeterhansen.org/class-notes/">here</a>. 