# In progress: Hansen, Khorrami & Tourre Example

````{prf:example}
This example is a discrete-time approximation of the model in Section 4.4 of {cite}`KhorramiTourreHansen:2024`. 
```{math}
C_t + I_t = \alpha K_t
```
Let  ${\widehat K}_t = \log K_t,$ ${\widehat C}_t = \log C_t$, and $\widehat{Z}_t^2 = \log Z_t^2$.  Then the endogenous state dynamics are given by:
```{math}
:label: equation1
\begin{align} 
{{\widehat K}_{t+\epsilon}} - {\widehat K}_t = \hspace{.2cm} & \left[ {\frac 1 \zeta}  \log \left( 1 + \zeta \frac{I_t}{K_t} \right) + \beta_k Z_t^1 - \eta_k \right]\epsilon  \cr \hspace{.2cm} &  - \frac{\epsilon}{2}|\sigma_k|^2 Z_t^2 \epsilon
 +  \sqrt{\epsilon Z_t^2}  {\sigma}_k  W_{t+\epsilon}  \end{align}, 
```
and the exogenous state dynamics by;
```{math}
:label: equation2
\begin{align}
Z_{t+\epsilon}^1 -
Z_{t}^1=   \hspace{.2cm} & - \beta_{1} Z_t^1 \epsilon +  \sqrt{\epsilon Z_t^2} \sigma_1  W_{t+\epsilon} \cr 
{\widehat Z}_{t+\epsilon}^2 - {\widehat Z}_{t}^2 =  \hspace{.2cm} &   - \beta_2 \left[ 1 - {\mu_2} \exp\left( - {\widehat Z}_t^2 \right) \right] \epsilon \cr \hspace{.2cm}&  - {\frac 1 2} |\sigma_2|^2  \exp\left( - {\widehat Z}_t^2 \right) \epsilon   
+ \sqrt{\epsilon} \exp\left( - {\frac 1 2}  {\widehat Z}_t^2 \right) \sigma_2 W_{t+\epsilon} \end{align} 
```
Observe that $W_{t+\epsilon}  $ is a multivariate standard normally distributed random vector.  

Let $X_t = [Z_t^1 \hspace{.2cm} {\widehat Z}_t^2],$ ${\widehat G}_t = {\widehat K}_t,  and $D_t = \frac{I_t}{K_t}.$ Rewrite resource constraint as:
```{math}
:label: equation3
{\widehat C}_t  = \log \left( \alpha - D_t \right) + {\widehat G}_t
```









The first-order conditions for the investment-capital ratio imply that  
```{math}
:label: equation4
\begin{align*}
0 &= - (1-\beta) \exp\left[(1- \rho) \left({\widehat C}_t - {\widehat G}_t \right)\right] \left( \frac{1}{\alpha - D_t} \right) 
\\&+ \beta \exp\left[(1 - \rho) \left({\widehat R}_t - {\widehat G}_t \right)\right] \left(\frac{1}{1 + \zeta D_t} \right)    
\end{align*}
```  
The co-state equation is:  
```{math}
:label: equation5
\begin{align*}
MX_t = \hspace{.2cm} & \beta \exp\left[(1 - \rho) \left({\widehat R}_t - {\widehat V}_t \right)\right] 
{\mathbb E} 
\left(
\exp\left[(1 - \gamma) \left({\widehat V}_{t+1} - {\widehat R}_t \right) \right] \psi_{x'}^x MX_{t+1} \mid {\mathfrak A}_t \right) 
\cr & \hspace{.2cm} + \beta \exp\left[(1 - \rho) \left({\widehat R}_t - {\widehat V}_t \right)\right] \psi_{x'}^g 
\end{align*}
```

```{math}
:label: equation6
\psi_{x'}^x = \begin{pmatrix} 1-\beta_1 \epsilon & \frac{1}{2} \exp\left(\frac 1 2 \widehat{Z}_t^2 \right) \sigma_1 \left(W_{t+\epsilon} - W_t \right) \\ 0 & 1 + \left(-\beta_2  \mu_2 + \frac 1 2 |\sigma_2|^2 \right) \exp(-\widehat{Z}_t^2) \epsilon - \frac{1}{2} \exp\left(- \frac 1 2 \widehat{Z}_t^2 \right) \sigma_2 \left(W_{t+\epsilon} - W_t \right) \end{pmatrix}
```
```{math}
:label: equation7
\psi_{x'}^g = \begin{pmatrix} \beta_k \epsilon & - \frac{1}{2}|\sigma_k|^2 \exp \left(\widehat{Z}_t^2\right) + \frac{1}{2} \exp\left(\frac 1 2 \widehat{Z}_t^2\right) \sigma_k\left(W_{t+\epsilon} - W_t \right)  \epsilon \end{pmatrix}
```