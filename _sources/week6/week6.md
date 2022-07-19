(week6)=
# Week 6: The connection between statistical mechanics and ensemble forecast I: Applications

## Connections between Liouville equation and ensemble forecasts
Last week, we learned the physical meaning of Jacobian matrix, determinant and Liouville equation. This week, we will walk through a few of their applications. Before we introduce the examples, we need to know there are two equivalent formulas of Liouville equation: (1) initial state-based Liouville equation and (2) current state-based Liouville equation. One easy way to distinguish these two is looking at how {eq}`eq44` is defined. 

First, let's consider a dynamical system only with first-order derivative. 
```{math}
:label: eq50
\begin{align}
\mathbf{\dot{X}} = \mathbf{\Phi}(\mathbf{X},t)
\end{align}
```
In the initial state-based method, the Liouville equation is only a function of initial state and time. Physically, we can know how a dynamical system evolves with time as long as we know its initial states and time for integration. i.e., 

```{math}
:label: eq51
\begin{align}
\mathbf{X} = \mathbf{X}(\mathbf{\Xi},t)
\end{align}
```
where $\Xi$ is the initial state and $t$ represent the number of time steps we integrate forward. Then the determinant of the Jacobian matrix, which exerts on the initial condition can be written as
```{math}
:label: eq52
\rho (\mathbf{\Xi},t) = \mathrm{det}(\frac{\partial \mathbf{X}}{\partial \mathbf{\Xi}})
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
(\frac{\partial }{\partial t} + \sum_{i}\Phi_i(\mathbf{X},t)\frac{\partial}{\partial \mathbf{X}})\rho = \rho(\mathrm{X},t)\frac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial \mathbf{X}}
```
compare {eq}`eq54` to the Liouville equation as a function of initial state (similar to {eq}`eq47`)

```{math}
:label: eq55
\frac{\partial \rho}{\partial t}  = \rho(\mathrm{\Xi},t)\frac{\partial \mathbf{\Phi}(\mathbf{X},t)}{\partial \mathbf{X}}\vert_{\mathbf{X(\mathbf{\Xi})}}
```
By observing {eq}`eq54` and {eq}`eq55`, it's not hard to find their analogs in equations of fluid dynamics. {eq}`eq54` is similar to so-called Eluerian form (i.e., observing the system at a fixed point), while {eq}`eq55` is similar to Lagrangian form (i.e, observers following the particle). This is also the reason why it's called "particle flow". 

## Liouville equation in 1st-order ODE case

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
\rho(\Xi,t) = \frac{\partial X}{\partial \Xi} = \frac{1}{(1+t\Xi)^{2}} 
```
and the Jacobian as a function of current state can be written as,  

```{math}
:label: eq60
\rho(X,t) = \frac{\partial X}{\partial \Xi} = \frac{1}{(1+t\frac{X}{1-X t})^{2}} = (1-tX)^2 
```













## References
```{bibliography} ../references.bib
:filter: docname in docnames
```

