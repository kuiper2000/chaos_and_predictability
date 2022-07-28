(week6)=
# Week 6: The connection between statistical mechanics and ensemble forecast I: Applications

## Connections between Liouville equation (LE) and ensemble forecasts
Last week, we learned the physical meaning of Jacobian matrix, determinant and LE. This week, we will walk through a few of their applications. Before we introduce the examples, we need to know there are two equivalent formulas of LE: (1) initial state-based LE and (2) current state-based LE. One easy way to distinguish these two is looking at how {eq}`eq44` is defined. 

First, let's consider a dynamical system with first-order derivative. 
```{math}
:label: eq50
\begin{align}
\mathbf{\dot{X}} = \mathbf{\Phi}(\mathbf{X},t)
\end{align}
```
In the initial state-based method, the LE is only a function of initial state and time. Physically, we can know how a dynamical system evolves with time as long as we know its initial states and time for integration. i.e., 

```{math}
:label: eq51
\begin{align}
\mathbf{X} = \mathbf{X}(\mathbf{\Xi},t)
\end{align}
```
where $\Xi$ is the initial state and $t$ represent the number of time steps we integrate forward. Then the determinant of the Jacobian matrix, which exerts on the initial condition can be written as
```{math}
:label: eq52
\alpha (\mathbf{\Xi},t) = \mathrm{det}(\frac{\partial \mathbf{X}}{\partial \mathbf{\Xi}})
```
While the dynamical system is deterministic, we can also find an 1-on-1 bijection between $\mathbf{\Xi}$ and $\mathbf{X}$ at a given t. i.e., 
```{math}
:label: eq53
\begin{align}
\mathbf{\Xi} = \mathbf{\Xi}(\mathbf{X},t)
\end{align}
```
The existence of 1-on-1 bijection enables us to rewrite {eq}`eq53` i.e., represent $\mathrm{\Xi}$ with $\mathrm{X}$ 
```{math}
:label: eq54
(\frac{\partial }{\partial t} + \sum_{i}\Phi_i(\mathbf{X},t)\frac{\partial}{\partial \mathbf{X}})\alpha = \alpha(\mathrm{X},t)\sum_{i}\frac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial \mathrm{X_i}}
```
compare {eq}`eq54` to the LE as a function of initial state (similar to {eq}`eq47`)

```{math}
:label: eq55
\frac{\partial \alpha}{\partial t}  = \alpha(\mathrm{\Xi},t)\sum_{i}\frac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial \mathrm{X_i}}\vert_{\mathbf{X(\mathbf{\Xi})}}
```
By observing {eq}`eq54` and {eq}`eq55`, it's not hard to find their analogs in equations of fluid dynamics. {eq}`eq54` is similar to so-called Eluerian form (i.e., observing the system at a fixed point), while {eq}`eq55` is similar to Lagrangian form (i.e, observers following the particle). This is also the reason why it's called "particle flow". 


To illustrative the mathematical equivalence between {eq}`eq54` and {eq}`eq55`, we use a Riccati equation (i.e., equation with form of $x^{'}=P(t)x^2+Q(t)x+R(t)$) to demonstrate 
```{math}
:label: eq56
\dot{X} = -X^2
```
the analytical solution is 
```{math}
:label: eq57
X = \frac{\Xi}{1+\Xi t}
```
where $\Xi$ is the initial condition and $t$ is time of integration. According to {eq}`eq52`, we can find the initial state as a function of current state. 

```{math}
:label: eq58
\Xi = \frac{X}{1-X t}
```
The Jacobian as a function of initial state can be written as, 

```{math}
:label: eq59
\alpha(\Xi,t) = \frac{\partial X}{\partial \Xi} = \frac{1}{(1+t\Xi)^{2}} 
```
and the Jacobian as a function of current state can be written as,  

```{math}
:label: eq60
\alpha(X,t) = \frac{\partial X}{\partial \Xi} = \frac{1}{(1+t\frac{X}{1-X t})^{2}} = (1-tX)^2 
```

Now, let's substitute $\alpha(\Xi,t)$ in {eq}`eq54` with {eq}`eq59` and $\alpha(X,t)$ in {eq}`eq55` with {eq}`eq60`, which will give us


```{math}
:label: eq61
\begin{align}
\frac{\partial \alpha}{\partial t}  = \alpha(\mathrm{\Xi},t)\frac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial \mathbf{X}}\vert_{\mathbf{X(\mathbf{\Xi})}} = \frac{1}{(1+t\Xi)^{2}} (-2) (\frac{\Xi}{1+\Xi t}) & = -2 X (\Xi,t) \alpha(\Xi,t) \\
(\frac{\partial }{\partial t} + \sum_{i}\Phi_i(\mathbf{X},t)\frac{\partial}{\partial \mathbf{X}})\alpha &= -2 X (1-tX)^2  
\end{align}
```
We can replace the right-hand-side of the second equation in {eq}`eq61` with {eq}`eq57`. Then we can find the first and second equations are mathematically equivalent. 


:::{note}
{eq}`eq51` leads to a very interesting debate in science, "Laplace demon theory". Laplace demon theory describes that if there is an intellect who knows the "specific location of every element particle in the universe", it will know the past as well as the future of the entire universe given {eq}`eq51` is deterministic. For example, every decision made by you and me can be considered as result of chemical reaction (i.e., the moving of electrons) in your brain. However, the discovery of uncertainty principle in 1927 protects us from having such demon.

In weather forecast, there is always some uncertainty since we can never have infinite small-scale observation.    
:::



## Solution of the Liouville equation 

Different from finding the solution of a dynamical system, the analytical solution of LE is relatively simple and interpretable. We can divide both side of {eq}`eq57` by $\alpha$ and integrate over $t$ dimension, which gives us, 


```{math}
:label: eq62
\int \frac{\partial \mathrm{ln}(\alpha)}{\partial t} dt = \int_{t'=0}^{t=\tau}\frac{\partial \mathbf{\Phi}(\mathbf{X},t')}{\partial \mathbf{X}}\vert_{\mathbf{X(\mathbf{\Xi},t')}} dt'
```
which further leads to 

```{math}
:label: eq63
\alpha(\mathbf{X},\tau) = \alpha(\mathbf{\Xi},0) \mathrm{exp}(\int_{t'=0}^{t=\tau}\frac{\partial \mathbf{\Phi}(\mathbf{X},t')}{\partial \mathbf{X}}\vert_{\mathbf{X(\mathbf{\Xi},t')}} dt')
```

or in density form 

```{math}
:label: eq64
\rho(\mathbf{X},\tau) = \rho(\mathbf{\Xi},0) \mathrm{exp}(-\int_{t'=0}^{t=\tau}\frac{\partial \mathbf{\Phi}(\mathbf{X},t')}{\partial \mathbf{X}}\vert_{\mathbf{X(\mathbf{\Xi},t')}} dt')
```

If we still remember what we leaned in [week1](https://kuiper2000.github.io/chaos_and_predictability/week1/week1.html#week1), {eq}`eq63` and {eq}`eq6` have exact the same formula. The only difference is, we leverage Taylor expansion to linearize the dynamical system while this week (and last week), we use Jacobian matrix to map the initial states to the final states. 

## Liouville equation in 1st-order ODE case

Here we use a simple ODE along with a python code showing how to implement LE in a forecast problem. 

```{math}
:label: eq65
\frac{d}{dt} x = x-x^3
```

Before we solve this ODE, we can take a look of its dynamical behavior. From {eq}`eq64`, one can easily find that it has three equilibrium point: -1, 0 ,1. At these three equilibrium points, $\frac{d}{dt}$ equals 0, indicating there is no acceleration. However, among these points, two of them are stable equilibrium points and the other one is unstable equilibrium point. 

The analytical solution of {eq}`eq64` can be written as (brainstorming...)
```{math}
:label: eq66
x(\Xi,t) = \Xi\mathrm{e}^{t} (1-\Xi^2+\Xi^2\mathrm{e}^{2t})^{-\frac{1}{2}}
```
where $\Xi$ is the initial state and $t$ is the time for integration. We can test if {eq}`eq56` is the solution of the differential equation {eq}`eq64` by substitute $t$ with $\infty$, which give us $x(\Xi,t) \rightarrow \pm 1$. This is consistent with the analytical analysis. 

On the other hand, {eq}`eq64` also leads to 

```{math}
:label: eq67
\frac{\partial \mathrm{\Phi}(x,t)}{\partial x} = 1-3x^{2}
```
.
by substituting {eq}`eq67` into {eq}`eq64`, we can have 
```{math}
:label: eq68
\rho(\mathbf{x},\tau) = \rho(\mathbf{\Xi},0) \mathrm{exp}(-\int_{t'=0}^{t=\tau} \{ 1-3(x \{ \Xi[x,t],t^{'} \})^{2} \}  dt')
```


Now, we can compare the probabilistic forecast generated by LE as well as the traditional forecast method. For the LE-based method, we can use numerical integration to approach the analytical solution: 

```python
def PDF(x,mu,std):
    return stats.norm(mu, std).pdf(x)

rho  = PDF(x,mu,std)

# x = location, rho = density in given location, n_step = number of time step, t_size = size of each step 
def LE_method(x,rho,n_step,t_size):
    for count in range(n_step): # forward for 1 time steps
        xi   = x
        xint = 0
        for i in np.arange(0,t_size*1000,1):
            tp = i*t_size/(t_size*1000)
            x  = xi*np.exp(tp)*(1-xi*xi+xi*xi*np.exp(2*tp))**(-1/2)
            c  = 1-3*x*x
            xint += c
        rho  = rho*np.exp(-xint*t_size/1000)
        rho  = rho/np.trapz(rho,x)
    return x,rho
```

For traditional method, we initialized model with 10000 normally distributed random points and integrate them forward in time. The comparison is shown below. 

```{figure} Liouville_vs_traditional.png
---
name: FIG13
---
The intercomparison between traditional probabilistic forecast method (bar) vs LE-based method.  
```

The LE-based method only used 200 data points (for the purpose of visualization) while the traditional method took 10000 scenarios to get similar answers, suggesting the LE-based method is 50 times more efficient than the traditional ensemble forecast method. 

## Liouville equation in Lorenz 63 model (Lorenz 84 model as HW)

We can also apply LE-based method to much higher dimension, such as Lorenz 84 model. Lorenz 84 model is a more complicated version of Lorenz 63 model (, which is also a simplification of 2D convection problem). 



```{math}
:label: eq68
\begin{align*}
\dot{X} &= -\sigma X+ \sigma Y \\
\dot{Y} &= -XZ+ rX- Y \\
\dot{Z} &= XY- bZ \\
\end{align*}
```  

or in a matrix form 

```{math}
:label: eq69
\mathbf{\dot{X}} = f(\mathbf{\dot{X}})
```  
To apply the LE-based method to the Lorenz 63 model, we first linearize the equation to simplify the problem given that the error growth rate is only the sum of its diagonal elements. (we will talk about the off-diagonal element in the chapters of generalized instability problem).  


```{math}
:label: eq70
\begin{bmatrix}
\dot{X} \\
\dot{Y} \\
\dot{Z} \\
\end{bmatrix} = 
\begin{bmatrix}
-\sigma & \sigma & 0 \\
-\bar{Z}+r & -1 & \bar{X} \\
\bar{Y} & \bar{X}  & -b
\end{bmatrix}
\begin{bmatrix}
X' \\
Y' \\
X' \\
\end{bmatrix}
```  

and the corresponding $\sum_{i}\frac{\partial \Phi_i(\mathbf{X},t)}{\partial \mathbf{X}}$ is

```{math}
:label: eq70
\sum_{i}\frac{\partial \Phi_i(\mathbf{X},t)}{\partial \mathbf{X}} = -\sigma -1 -b
``` 

Interestingly, we can find the error growth (i.e., decaying of $\rho$) is independent of location and this seems counter-intuitive? Lorenz 63 system is well-designed model, where each component has a constant error growth rate. However, it doesn't mean the PDF will stay Gaussian given that each member are moving around within a nonlinear system. Physically, ensemble members are something like particles carrying different amount of electric charges (probability density). Thus, even each member has similar shape of PDF (Gaussian), their superimposition in phase space are not necessary Gaussian.  

The animation below shows the LE-based method of probabilistic forecast for Lorenz 63 model. We can find that the system goes through a stage of linear error growth first (stretching along y-axis). At the later stage, we can find the regions with higher probability are concentrated around two attractors. This feature is generally consistent with our finding. 


```{figure} Lorenz_LE.gif
---
name: FIG14
---
Using LE-based method of generating probabilistic forecast in Lorenz 63 model. 
```


:::{note}
In the homework assignment, you will practice how to apply the LE-based method to Lorenz 84 model ({eq}`eq71`). Different from the Lorenz 63 model, the convergence/divergence of $\Phi_i(\mathbf{X},t)$ is state-dependent (a function of location). This indicates that the probability density carried by each member might experience nonlinear decaying/growth when it is moving around the dynamical system.  

```{math}
:label: eq71
\begin{align*}
\dot{X} &= -Y^{2}-Z^{2}-aX+aF \\
\dot{Y} &= XY-bXZ-Y+G \\
\dot{Z} &= bXY+XZ-Z \\
\end{align*}
```
:::









## References
```{bibliography} ../references.bib
:filter: docname in docnames
```

