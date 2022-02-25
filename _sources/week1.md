# Week 1: Predictability of Weather and Climate

In this section, we will walk through a few fundamental concepts of predictability, including how does it connect to statistics, differential equations and linear algebra. You will find out the definition of predictability is relatively straightforward. 

__Climatological distribution and forecast distribution__
\
When we say something (let's say $x(t)$ in our case) is predictable, that means we can somewhat track its time evolution. However, in most cases, our confidence about what $x(t)$ should be will decrease with the increase of forecast lead time (defined as the time difference between initial state and the final state). Thus, it is physically more reasonable to use a probability density function (PDF $p(x(t))$) to describe $x(t)$ since we are talking about how "confident" we are. 
\
\
Now, let's imagine two very extreme cases. The first case is $t\rightarrow 0$. In this case, we are very certain what $x(t)$ should be because it is literally the current states or so-called _nowcasting_. Thus, it's not hard to find $p(x(t))$ is very narrow and almost like a delta function ([FIG1](FIG1), solid line). The second case is that $t\rightarrow \infty$. In this case, we don't have any confidence. Thus, the best way to estimate $x(t)$ is by randomly sampling from its historical values (Fig. 1, dashed line, a.k.a. guessing). Thus, from $t=0$ to $t\rightarrow \infty$, we can find $p(x(t))$ is evolving from one PDF to the other. Here we give these two distribution names. For the one which we use to describe how confident we are is called _forecast distribution_ and the one when forecast lead time approaches infinity is called _climatological distribution_ (it's also called "population" in statistics). 

```{figure} week1/FIG1.png
---
name: FIG1
---
An example of forecast probability density function and climatological probability density function. 
```
\
Here we can formulte our mathematical definition of "time of predictability limit". It is defined as the "moment" that the following null hypothesis is rejected.

```{math}
:label: eq1

H_o: p(x(t))\neq p(x(\infty)) 
```
\
{eq}`eq1` indicates that the moment we can no longer tell the difference between $p(x(t))$ and $p(x(\infty))$ is also the moment we hit the predictability limit because the best estimation of $p(x(t))$ is not better than random guess! Now, we can see how the idea of predictability connects to statistical test. 
\
Another interesting thing you might have noticed...when we talked about predictability, we didn't use the ground truth, i.e., the observed $x(t)$. Yes, because the measurement of predictability doesn't rely on ground truth or observation. Instead, it only relis on the forecast states. This is so-called _perfect model assumption_ and we will have more detailed discussion in week 3.   
\
__State-dependent predictability and mathematical assumptions__
\
When talking about predictability, one should notice there is no universal value for predictability limit. You might think...OK that sounds weird and counterintuitive. For example, we know for typical numerical weather forecast, we will say the predictability limit is around 10 days to 2 weeks. After that, we can no longer trust the model output. However, you will soon realize "10 days to 2 weeks" is just rule of thumb to describe the averaged predictability in current NWP systems. In some cases, we can even struggle with the low prediction confidence at a forecast lead time of 3 days! Sounds crazy right? Here, let me use the famous Lorenz 63 model to demonstrate what I mean. 

:::{note}
Lorenz 63 model can be considered as the minimalist model for studying chaos and predictability and it only contains 3 progonostic variables. Indeed, to generate chaotic behavior, the minimum size of independent variables is 3. Smaller than this number, we won't have chaotic behavior in a dynamical system. We will cover more details about Lorenz model in week 3. 
:::

```{math}
:label: eq2
\begin{align*}
\frac{dx}{dt} &= \sigma (y-x) \\
\frac{dy}{dt} &= x(\rho-z)-y  \\
\frac{dz}{dt} &= xy-\beta z 
\end{align*}
```
\
{eq}`eq2` is the Lorenz 69 model with the selected parameters of $\rho=28$, $\sigma=10$ and $\beta=\frac{8}{3}$. This parameter selection makes sure Lorenz 69 model has a chaotic behavior. (The homework assignment will let you test different parameter selections and see how the dynamical behavior changes.) In this model, we select three different initial states and generate the ensemble simulations by slightly perturbing the initial x,y and z. The result is shown in Fig. 2. We can find in [FIG2(a)](FIG2) case, the model spread grows relatively fast and soon diverge into two groups ($\sim 50\%$ on the right/left). In the second scenario (i.e., [FIG2(b)](FIG2)), the difference among ensemble members is small at the very beginning but start to bifurcate in the second-half of our simulations. In the last case ([FIG2(c)](FIG2)), different members stay close even at the end of simulation perod (red dots) suggesting a relatively small growth rate.   

```{figure} week1/FIG2.png
---
name: FIG2
---
The three scenarios of ensemble forecasts based on L69 model (a) fast error growth (b) average error growth (c) slow error growth
```
We can go one step further and plot the forecast PDF at the final state from three scenarios and the result is shown in [FIG3](FIG3). We also include the climatological PDF as a reference, which is generated by a long-term simulation with 50000 time steps (non-dimensional time). If you still remember what we discussed at the very beginning, the predictability limit is defined by "when" we can reject the null-hypothesis in {eq}`eq1`. Thus, by comparing how similar the two distributions are can inform us how predictable the given dynamical system is. In [FIG3](FIG3), we can easily find the forecast PDF in (a) and the climatological PDF are quite similar. However, for cases in (b) and (c), we can still distinguish the forecast PDFs from the climatological PDF indicating we haven't reached the predictability limit yet.  

```{figure} week1/FIG3.png
---
name: FIG3
---
The final states' probability density function of three ensemble forecast shown in Fig. 2
```
One take-home message in [FIG2](FIG2) and [FIG3](FIG3) is, the predictability is a function of state rather than an universal number. Now, we can imagine why the most experienced forecasters and the state-of-art NWP systems can still struggle sometimes...(Thus, be kind to them). There are also a few things to keep in mind. First, we only consider the uncertainty in model initial states but have ignored the uncertainty in model structures (e.g., if the model is correct) or the rounding error. This _perfect model assumption_ is one of the most important assumption in the enrire course (and probably the entire field) since it tells us the _upper limit_ of predictability. It also tells us, as long as we have infinitesimal error in model initial states, the predictability limit is inevitable. Second, _predictability_ has no meaning unless the ensemble simulations are used since its definition is based on how fast one ensemble member can diverge from the other. Third, to reject the null hypothesis in {eq}`eq1`, we need to choose a significance level. This indicates that we can always use a very low significance level (e.g., $1\%$ or $0.1\%$) to say we haven't reached the predictability limit. Thus, clearly define what threshold we are using for testing the null hypothesis is very important. 

__Where the uncertainties come from?__ 
\
At the end of previous section, we talk about the potential uncertainties can come from three different places (or at least we can attribute any kind of forecast uncertainty to these three). Of course they can happen at the same time and be indistinguishable in some cases but we will first scrutinize them individually for the purpose of discussion. The first uncertainty is the initial state error or the observational error. In a perfect observation, the spatial and temporal resolution should go all the way down to the smallest scale (molecular scales). Bcause as long as we have missing observation, the upscale growth of initial error from those regions will ultimately lead to an unpredictable future (if the underlying dynamics is chaotic). 
\
\
The second error source is from the imperfect model physics. Specifically, the physical parameterizations used to approximate the bulk effect of subgrid-scale processes (i.e., the scales smaller than the model grid) are the major uncertainty source. The details of parameterizations will be discussed in other class and we will briefly walk through the main concept. One reason of using physical parameterizations is due to the limited computational power. For example, to accucrately predict the time evolution of an extratorpical storm, we also need to predict the convection imbeded in the frontal structures since the latent heat release by these convections is not negligible. However, explicitly resolving this small-scale thunderstorms is not computationally feasible for the purpose of synoptic weather forecast. Thus, in most NWP systems, we use so-called cumulus parameterizations to approximate the bulk effect of convective cloud. The reason that cumulus parameterization works is that the large-scale environment usually reach an quasi-equilibrium state with the small scale convection (i.e., coherence exists). Therefore, we can approximate the net convective activity by using the large-scale information. However, similar to the first uncertainty (observational uncertainty), the observational error exists in all scales and thus a PDF (stochastic parameterization) is required to describe the subgrid-scale statistics as well rather than a single value. While it seems necessary, it's not the case in most prevailing NWP systems. In addition to the cumulus parameterization, similar problems exist in boundary layer parameterizations.      
\
The last uncertainty is the rounding error. Comparing with the former two uncertainties, rounding error has the least impacts to the weather and climate predictions. Sometimes we can even have some trade-off... i.e., allowing for certain rounding error to save some computational time. The main reason we can do that is because the uncertainties from observation and model physics are way bigger ($>\mathcal{O}(5)$) than the uncertainty of rounding error. More details can be found {cite}`hatfield2020single`.

:::{note}
Lorenz 96 is one of the simplest models attemping to deal with the underpinning theory of subgrid-scale processes, and we will talk about more details in Week 4.  
:::
\
__Introduction to ensemble forecast in weather and climate__ 
Given the discussion about, we know (1) ensemble forecast is necessary when talking about predictability (2) predictability limit is state-dependent. We also use Lorenz 63 model to reaffirm these two conclusions. Now, let's move to a more generalized form and try to include both conclusions in  

```{bibliography}
:filter: references
```
