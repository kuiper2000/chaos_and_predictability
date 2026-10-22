(week7)=
# Week 7: The stability theory

In week 6, we learn how to estimate the time evolution of PDF. We can easily find that the analytical solution of LE (i.e., {eq}`eq64`) follows a form of exponential function, where the integral inside the exponent represents the state-dependent error growth rate. For example, if the local divergence $\sum_{i}\frac{\partial \Phi_i}{\partial X_i}$ is positive, the exponent is negative, indicating the volume (density) spanned by ensemble members will increase (decrease) with time at that given location. Thus, it's pretty much like a stability problem we usually deal with at the traditional geophysical fluid dynamics. For example, in mid-latitude, some regions with strong meridional temperature gradient (vertical wind shear), we can see frontal genesis. These regions have very strong baroclinicity and thus a small disturbance can easily grow. In atmospheric dynamics, we can calculate so-called Eady growth rate, ({eq}`eq74`, where $f_0$ is the Coriolis parameter, $N$ is the Brunt-Vaisala frequency (i.e., the restoring force by buoyancy), $H$ is the depth of the layer, and $\Delta U$ is the change of zonal wind across that layer, so that $\frac{\Delta U}{H}$ is the vertical wind shear) to estimate the environment favoring the baroclinicity.    

```{math}
:label: eq74
\sigma = 0.31\frac{f_0}{NH}\Delta U
``` 
{eq}`eq74` is a convenient formula for estimating if an environment is baroclinically unstable. However, there are many different kinds of instability and thus a more general framework is required to quantify if an error in initial states will grow with time. In this week, we will introduce three traditional methods (1) Direct method, (2) Energy method, and (3) von Neumann method. All of these methods are widely used in numerical and dynamical system analyses. If any of them indicates instability, the system is either numerically or physically unstable. (One should notice that the system can be physically stable but numerically unstable). Here we use a very simple 1D-advection model to demonstrate how to apply three different methods to analyze the stability in a given system. 


```{figure} 1D_advection.png
---
name: FIG15
---
An example of a one-dimensional grid, with uniform grid spacing $\Delta x$. The grid points are denoted by the integer index j. Half-integer grid points can also be defined, as shown by the x at j+$\frac{1}{2}$. 
```

We first start with grid spacing. [FIG15](FIG15) is schematic diagram of 1D grid (sometimes we use mesh or lattice). Since each grid is equally spaced, we can define a function $f(x_j)=f(j\Delta x)$, which further leads to three different definitions of finite difference: (1) For forward difference at point $j$, we have $\frac{f(x_{j+1})-f(x_{j})}{\Delta x}$, (2) for backward difference at point $j$, we have $\frac{f(x_{j})-f(x_{j-1})}{\Delta x}$ and (3) for central difference at point $j$, we have $\frac{f(x_{j+1})-f(x_{j-1})}{2\Delta x}$ (note the $2\Delta x$ in the denominator, because the two points straddling $j$ are $2\Delta x$ apart). 


Now, we can take a look of the 1D-advection model. 

```{math}
:label: eq75
\frac{\partial A}{\partial t} + c (\frac{\partial A}{\partial x})= 0   
``` 

with an initial condition

```{math}
:label: eq76
A(x,0) = F(x) 
``` 

Physically, {eq}`eq75` carries $A_{j}^{n}$ at a position of a particle with a constant speed $c$. We interpret $A_{j}^{n}$ as a conserved property of the particle. To find the analytical solution, we can assume $\xi \equiv x-ct$. A line of constant $\xi$ is then a particle trajectory: holding $\xi$ fixed and differentiating gives $dx-cdt=0$, i.e., $(\frac{\partial x}{\partial t})_{\xi}=c$ and equivalently $(\frac{\partial t}{\partial x})_{\xi}=\frac{1}{c}$. Using the chain rule, we further find that 

```{math}
:label: eq77
\begin{align*}
(\frac{\partial A}{\partial x})_{\xi} 
&= (\frac{\partial A}{\partial x})_t(\frac{\partial x}{\partial x})_{\xi}+(\frac{\partial A}{\partial t})_x(\frac{\partial t}{\partial x})_{\xi}  \\
&= (\frac{\partial A}{\partial x})_t+(\frac{\partial A}{\partial t})_x(\frac{1}{c})  \\ 
&= \frac{1}{c}\left[(\frac{\partial A}{\partial t})_x+c(\frac{\partial A}{\partial x})_t\right] \\
& = 0
\end{align*}
``` 

and 

```{math}
:label: eq78
\begin{align*}
(\frac{\partial A}{\partial t})_{\xi} 
&= (\frac{\partial A}{\partial t})_x(\frac{\partial t}{\partial t})_{\xi}+(\frac{\partial A}{\partial x})_t(\frac{\partial x}{\partial t})_{\xi}  \\ 
&= (\frac{\partial A}{\partial t})_x+(\frac{\partial A}{\partial x})_t c \\
& = 0
\end{align*}
``` 

Note that everything up to the second line of {eq}`eq77` and {eq}`eq78` is pure calculus; it holds for any field $A(x,t)$. The physics only enters in the last line, where the bracketed expression is recognized as the left hand side of {eq}`eq75` and is therefore zero. Note also that the two results are not independent: because $dx=cdt$ along a characteristic, {eq}`eq77` is simply {eq}`eq78` divided by $c$ (assuming $c \neq 0$).

{eq}`eq77` and {eq}`eq78` suggest that $A$ will conserve as long as it stay on the characteristic line $\xi$. Thus, $\xi$ is an analytical solution of {eq}`eq75` (i.e., $A=f(\xi)$). Now let's consider using a numerical way to solve {eq}`eq75`, which can be written as the form below,


```{math}
:label: eq79
\frac{A^{n+1}_{j}-A^{n}_{j}}{\Delta t} + c (\frac{A^{n}_{j}-A^{n}_{j-1}}{\Delta x})= 0   
``` 

Here we use forward difference quotient for time and backward difference for space. We can easily find that what determines the time tendency is the information from _upstream_ (i.e., information is only from one side, i.e., $A_{j-1}^{n}$, and has no dependence on the other side, i.e., $A_{j+1}^{n}$). Let $A(x_j,t_n)$ be the exact solution at the point $x=j\Delta x$ and $t=n\Delta t$, and let $A^{n}_{j}$ be the numerically solved solution at that same grid point. In general, $A(x_j,t_n) \neq A^{n}_{j}$ in most cases (but we definitely hope they are the same ;) ). Next, we will examine if {eq}`eq79` is numerically stable. 

First, we introduce something called _truncation error_. If we plug the exact solution $A(x_j,t_n)$ into the continuous equation {eq}`eq75`, the left hand side is 0 by definition. However, if we plug that same exact solution into the _discrete_ operator of {eq}`eq79`, the result won't be 0 for most cases, i.e., 

```{math}
:label: eq80
\frac{A(x_{j},t_{n+1})-A(x_{j},t_{n})}{\Delta t} + c (\frac{A(x_{j},t_{n})-A(x_{j-1},t_{n})}{\Delta x})= \epsilon   
``` 
and applying Taylor expansion, we know. 

```{math}
:label: eq81
\epsilon = (\frac{1}{2!}\Delta t \frac{\partial^2 A}{\partial t^2}+\ldots)+c(-\frac{1}{2!}\Delta x\frac{\partial ^2 A}{\partial x^2}+\ldots)  
``` 

We say this is a first order scheme given that the leading terms of $\epsilon$ are proportional to the first power of $\Delta t$ and $\Delta x$. To reduce the _truncation error_, we can select both $\Delta t$ and $\Delta x$ as small as possible. 

## Direct Method

In addition to the truncation error which estimates how accurate a numerical differential equation is, we also have _discretization error_, which estimates how accurate the answer of a differential equation is. However, a decrease in the truncation error does not necessarily ensure a decrease in the discretization error. Discretization error is the difference between the exact solution, $A(x_j,t_n)$, and the discretized solution $A_{j}^{n}$. An easy way to visualize the difference between analytical solution and the discretized solutions is using so-called _domain of dependence_. i.e., plot the discretized solution as a function of $n$ and $j$. 


```{figure} domain_dependence.png
---
name: FIG16
---
A schematic diagram of domain dependence. The shading region shows the "domain of dependence" of {eq}`eq79`, i.e., the set of earlier grid points that $A^{n}_{j}$ actually depends on.  
```

The shading regions in [FIG16](FIG16) show where the numerical solution should fall in order to get smaller discretization error. If today, the integration starts at a point of $A(x_0,0)$ where $x_0<0$ as indicated by the thin line in [FIG16](FIG16). We can find that the discretization error will never decrease as long as $\frac{\Delta x}{\Delta t}$ is a constant. 

To understand what happens here, we can leverage the analytical solution of {eq}`eq75`. For the numerical simulation starts from $A(x_0,0)$, its characteristic line is $\Delta x-c\Delta t=x_0$. Rearranging the equation, it gives us $(1-\frac{x_0}{\Delta x})=c\frac{\Delta t}{\Delta x}$. 

```{math}
:label: eq82
(1-\frac{x_0}{\Delta x})=c\frac{\Delta t}{\Delta x}
``` 

We can easily find, if the discretized solution is identical to the analytical solution, then the second term on the left hand side of {eq}`eq82` will be 0. Therefore, the second term on the left hand side can be considered as the _discretization error_, which measure the difference between analytical and numerical solutions. In other words, when $c\frac{\Delta t}{\Delta x}=1$, we are certain that $x_0$ is 0. Obviously, $c\frac{\Delta t}{\Delta x}$ is an important parameter for quantifying how close the two solutions (numerical and analytical) are. This parameter is called the Courant number, and the stability constraint built on it is the CFL criterion, named after three great mathematicians: Courant, Friedrichs and Lewy. 

An alternative way to interpret {eq}`eq82` is through "interpolation and extrapolation". For example, we can rewrite the upstream scheme as 

```{math}
:label: eq83
A^{n+1}_{j} = (1-\mu) A^{n}_{j} + \mu A^{n}_{j-1}
``` 

where $\mu$ equals $c\frac{\Delta t}{\Delta x}$. From the schematic diagram below, we can know the solution in the next time step will always be bounded within $A_{j}^{n}$ and $A_{j-1}^{n}$ if $\mu$ is smaller than 1 and bigger than 0. Thus, even in the worst scenario, the signals will be smoothed out rather than exploding. 

```{figure} interprolation_extrapolation.png
---
name: FIG17
---
A schematic diagram of "interpolation vs extrapolation"
```

We can rewrite {eq}`eq83` in an absolute form (by the triangle inequality) and the following conclusion holds as long as $0\leq\mu \leq 1$, so that both coefficients $1-\mu$ and $\mu$ are non-negative

```{math}
:label: eq84
|A^{n+1}_{j}| \leq (1-\mu) |A^{n}_{j}| + \mu |A^{n}_{j-1}|
``` 

and

```{math}
:label: eq85
\begin{align*}
\mathrm{max}_{(j)} |A^{n+1}_{j}| &\leq (1-\mu) \mathrm{max}_{(j)}|A^{n}_{j}| + \mu \mathrm{max}_{(j)}|A^{n}_{j-1}| \\
&= \mathrm{max}_{(j)}|A^{n}_{j}|
\end{align*}
```  

The second line of {eq}`eq85` follows because of $\mathrm{max}|A_{j}^{n}| = \mathrm{max}|A_{j-1}^{n}|$ (shifting the index does not change the maximum over all $j$), so that $(1-\mu)+\mu=1$. In other words, the largest amplitude in the domain can never grow from one time step to the next. By inspecting {eq}`eq85`, we can again find that the criteria for a numerically stable scheme is 


```{math}
:label: eq86
0\leq \mu \leq 1
```  

This approach gives us a direct measurement of the instability in a numerical scheme or any phenomenon of interest. It's called _direct method_. 

## Energy Method

The second approach is so-called _energy method_. Mathematically, we can repeat the same exercise to acquire the same conclusion. I will leave it to the readers to derive the formula. Details can be found in homework instruction. 

## von Neumann Method

In this part, we will walk through another very powerful tool for testing instability, which is called von Neumann method (named after a Hungarian mathematician, John von Neumann). Let's return to {eq}`eq75` and assume the domain size is infinite for simplicity. This will lead to a wave solution  

```{math}
:label: eq87
A(x,t) = \mathrm{Re}[\hat{A}(t)e^{ikx}]
```  

where $|\hat{A}|$ is the wave amplitude. Here we consider a simple wave number for convenience but we can replace the right hand side of {eq}`eq87` with a summation over different wave numbers. Substituting {eq}`eq87` into {eq}`eq75`, we have


```{math}
:label: eq88
\frac{d \hat{A}(t)}{dt} + ikc \hat{A}(t) = 0
```

and the corresponding solution is 

```{math}
:label: eq89
\hat{A}(t) = \hat{A}(0)e^{-ikct}  
```
where $\hat{A}(0)$ is the initial states. We can substitute this "time-only" part back to {eq}`eq87` to derive the whole solution (i.e., with both space and time structure). 

```{math}
:label: eq90
A(x,t) = \mathrm{Re}[\hat{A}(0)e^{ik(x-ct)}]
```

For a finite difference problem, we can assume that the solution has a form of 

```{math}
:label: eq91
A^{n}_{j} = \mathrm{Re}[\hat{A}^{n}e^{ikj\Delta x}]
```
where $\hat{A}^{n}$ is the amplitude of wave $e^{ikj\Delta x}$ at the nth time step. One should notice that $\hat{A}^{n}$ can be an imaginary number, which determines both the phase angle and amplitude of a wave at a given moment. 

Now we can define an amplification factor, which measures the growth of $A$. 
```{math}
:label: eq92
\hat{A}^{n+1} \equiv \lambda \hat{A}^{n}
```

For the analytical solution of {eq}`eq75`, advancing {eq}`eq89` by one time step gives 

```{math}
:label: eq93
\hat{A}(t+\Delta t) = e^{-ikc\Delta t}\hat{A}(t) 
```

so the amplification factor is $\lambda = e^{-ikc\Delta t}$. This implies that for the exact advection problem with a single wave number, $|\lambda| = 1$, i.e., the signal won't amplify or decay; $\lambda$ only rotates the phase. 



The question is: _what is amplification factor in discrete form?_ (i.e., numerical solution). If we start from $t=0$, the numerical solution can be written as a form of amplification factor.

```{math}
:label: eq94
\hat{A}^{n} \equiv \hat{A}^{0} \lambda^{n}
```

therefore, the requirement for model not to blow up is $|\lambda|\leq 1$. To explicitly derive the discrete form of amplification factor, we can substitute {eq}`eq91` back to {eq}`eq79`. Every term then carries a common factor $e^{ikj\Delta x}$, which we divide out (note that $A^{n}_{j-1}$ contributes $e^{ik(j-1)\Delta x}=e^{ikj\Delta x}e^{-ik\Delta x}$), and have


```{math}
:label: eq95
\frac{\hat{A}^{n+1}-\hat{A}^{n}}{\Delta t}+(\frac{1-e^{-ik\Delta x}}{\Delta x})c\hat{A}^{n} = 0 
```

we can find that $\hat{A}^{n+1}=\left[1-c\frac{\Delta t}{\Delta x}(1-e^{-ik\Delta x})\right]\hat{A}^{n}$. Comparing with the {eq}`eq88`, it is found $\frac{1-e^{-ik\Delta x}}{\Delta x}c$ taking place of $ikc$. To prove this point, we can force $\Delta x \rightarrow 0$ and therefore $\frac{1-e^{-ik\Delta x}}{\Delta x}c \sim ikc$. 


If we use the definition of {eq}`eq92` and {eq}`eq95`, we can find 

```{math}
:label: eq96
\lambda = 1-\mu (1-\mathrm{cos} k\Delta x+i\mathrm{sin} k\Delta x)
```

where $\mu$ is the Courant number and we used $e^{-ik\Delta x}=\mathrm{cos}\,k\Delta x-i\,\mathrm{sin}\,k\Delta x$. {eq}`eq96` is what we are looking for. However, it's in a complex form. Thus, we multiply {eq}`eq96` by its complex conjugate to derive the squared modulus of {eq}`eq96`. 

```{math}
:label: eq97
|\lambda|^{2} = 1+2\mu(\mu-1)(1-\mathrm{cos}(k\Delta x))
```

Since $1-\mathrm{cos}(k\Delta x)\geq 0$ for every wave number, the requirement $|\lambda|^{2}\leq 1$ reduces to $\mu(\mu-1)\leq 0$, i.e., $0\leq \mu \leq 1$. This is exactly {eq}`eq86`, which we obtained from the direct method: two very different routes, the same CFL criterion.

Obviously, the amplification factor is also a function of wave number. Physically this make senses because when we say that a signal grows, not only its amplitude but also its phase angle should be considered. For example, when a vector rotates from $[1,0]$ to $[\frac{\sqrt{2}}{2},\frac{\sqrt{2}}{2}]$, the y component increases but the overall amplitude (i.e., $\sqrt{x^{2}+y^{2}}$) doesn't change. Different wave numbers correspond to different rotation rate. We will cover more details in the next section when we talk about the nonmodal system. 

## Stability in Geophysical Fluid Dynamics 

The stability theory has been widely applied to different fields. What we learn above is one application in numerical modeling. Here, we leverage barotropic instability for demonstrating how to apply those tools to GFD problems. 



Let's setup our question first. Imagine there is a 2D barotropic flow, which is bounded in meridional direction and periodic in zonal direction. Around the center of domain, we prescribe a background flow (i.e., [FIG18](FIG18)), which is positive $\bar{U}$ when $y>L$, is negative $\bar{U}$ when $y<-L$, and varies linearly as $\bar{U}\frac{y}{L}$ in between. In addition, no zonal structure of background flow is given for simplification. Note that this profile is continuous everywhere but its shear $\bar{U}_y$ jumps at $y=\pm L$, which is exactly where the interesting physics will come from. 

```{figure} Shear_instability.png
---
name: FIG18
---
A schematic diagram of prescribed background flow for shear instability. 
```

The corresponding governing equations can be written as:

```{math}
:label: eq98
\begin{align}
\frac{\partial u}{\partial t} +\bar{U}\frac{\partial u}{\partial x}+v\frac{\partial \bar{U}}{\partial y} &= -\frac{1}{\rho_0}\frac{\partial p}{\partial x} \\
\frac{\partial v}{\partial t} + \bar{U}\frac{\partial v}{\partial x} &=-\frac{1}{\rho_0}\frac{\partial p}{\partial y} \\
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y} &= 0 
\end{align}

``` 

Consider solutions of the form 

```{math}
:label: eq99
\begin{align}
u &= \tilde{u}(y)e^{ik(x-ct)} \\
v &= \tilde{v}(y)e^{ik(x-ct)} \\
p &= \tilde{p}(y)e^{ik(x-ct)}
\end{align}
``` 

where c is a complex, so that a positive imaginary part $c_i>0$ means the wave grows like $e^{kc_i t}$. Substitute into {eq}`eq98` gives (with $\alpha_0 \equiv \frac{1}{\rho_0}$)

```{math}
:label: eq100
\begin{align}
ik(\bar{U}-c)\tilde{u}+\bar{U}_y\tilde{v} &= -\alpha_0 ik\tilde{p} \\
ik(\bar{U}-c)\tilde{v} &= -\alpha_0 \tilde{p}_y \\
ik\tilde{u} + \tilde{v}_y & = 0
\end{align}
``` 

by cross-differentiating, we can eliminate $\tilde{p}$ and $\tilde{u}$ and have a single O.D.E. in $\tilde{v}$ 

```{math}
:label: eq101
\tilde{v}_{yy}-\tilde{v}(k^2+\frac{\bar{U}_{yy}}{\bar{U}-c}) = 0
``` 

Away from $y=\pm L$ the background flow is linear in $y$, so $\bar{U}_{yy}=0$ and {eq}`eq101` reduces to $\tilde{v}_{yy}-\alpha^{2}\tilde{v}=0$ with $\alpha \equiv |k|$. The general solution for y structure of meridional wind is therefore 

```{math}
:label: eq102
\begin{align}
\tilde{v} & = Ae^{-\alpha y}            &\; \mathrm{if} \; y\geq L  \\
\tilde{v} & = Be^{-\alpha y}+Ce^{\alpha y}   &\; \mathrm{if} \; L\geq y \geq -L \\
\tilde{v} & = Fe^{\alpha y}             &\; \mathrm{if} \; -L\geq y 
\end{align}
```   

The growing exponential has already been discarded in the two outer regions so that $\tilde{v}\rightarrow 0$ as $y\rightarrow \pm \infty$.


Now, we need some boundary conditions to solve $A,B,C$ and $F$ in $\tilde{v}$. We have already imposed $\tilde{v}=0$ when $y\rightarrow \pm \infty$. We also need some other conditions to solve the O.D.E. given that the background flow is not differentiable at $y=\pm L$. The first one is the continuity of pressure (a jump in pressure would imply an infinite acceleration of the fluid there). Eliminating $\tilde{u}$ between the first and third rows of {eq}`eq100` gives $ik\alpha_0\tilde{p}=(\bar{U}-c)\tilde{v}_y-\bar{U}_y\tilde{v}$, so continuity of $\tilde{p}$ gives us 

```{math}
:label: eq103
[(\bar{U}-c)\frac{d \tilde{v}}{dy}-\bar{U}_y\tilde{v}]^{y^{+}_{0}}_{y^{-}_{0}} = 0
```   

Also, the displace of fluid should be continuous. i.e., $\tilde{v}\equiv\frac{d}{dt}\delta y$ where $\delta y$ is the displacement. 

```{math}
:label: eq104
\begin{align}
\tilde{v} & \equiv\frac{d}{dt}\delta y  \\
& = (\frac{\partial }{\partial t}+\bar{U}\frac{\partial}{\partial x})\delta y \\
& = (-ikc+ik\bar{U})\delta y
\end{align}
``` 
because $\bar{U}$ is continuous over the entire domain and $\delta y$ is continuous, that implies $\tilde{v}$ is also continuous. Thus, {eq}`eq104` will reduce to 

```{math}
:label: eq105
Ae^{-\alpha L}=Be^{-\alpha L}+Ce^{\alpha L} \\
Fe^{-\alpha L}=Be^{\alpha L}+Ce^{-\alpha L}
``` 

Also, apply the jump condition in {eq}`eq103` at $y=L$ and $y=-L$, we can have another two equations (recall $\bar{U}_y=0$ in the two outer regions and $\bar{U}_y=\frac{\bar{U}}{L}$ inside the shear layer)
```{math}
:label: eq106
-\alpha A(\bar{U}-c)e^{-\alpha L}=[-\alpha B(\bar{U}-c)-B\frac{\partial \bar{U}}{\partial y}]e^{-\alpha L}+[\alpha C(\bar{U}-c)-C\frac{\partial \bar{U}}{\partial y}]e^{\alpha L} \\
\alpha F(-\bar{U}-c)e^{-\alpha L}=[-\alpha B(-\bar{U}-c)-B\frac{\partial \bar{U}}{\partial y}]e^{\alpha L}+[\alpha C(-\bar{U}-c)-C\frac{\partial \bar{U}}{\partial y}] e^{-\alpha L} 
``` 
Four equations with four unknowns. For a nontrivial solution (i.e., $A,B,C,F$ not all zero), the determinant of the system must vanish, and we can find the following dispersion relationship. 

```{math}
:label: eq107
(\frac{2\alpha L c}{\bar{U}})^2= (2\alpha L-1)^{2}-e^{-4\alpha L}
``` 

Ok, after so many steps of derivation, you might wonder where is the stability analysis. If you still remember von Neumann method, we first assume a wave solution and then examine if the solution will grow or not to determine the key parameters. That's exactly what we are doing for this shear stability analysis! Here the key parameter is $\alpha L$, the wave number scaled by the width of the shear zone. When the right hand side of {eq}`eq107` is negative, $c$ is purely imaginary, one of the two roots has $c_i>0$, and the wave grows: the flow is unstable. That happens for $\alpha L < 0.6392$. In other words, only waves that are _long compared with the width of the shear zone_ can extract energy from the background flow; short waves do not feel the two edges at $y=\pm L$ at the same time and are neutral. This short-wave cutoff is the barotropic analogue of the CFL criterion we found for the numerical scheme: a single dimensionless number decides growth or decay.

## Howard's Semicircle Theorem

{eq}`eq107` is the exact answer, but it took four boundary conditions and a $4\times4$ determinant to get, and it only applies to this one background flow. It is natural to ask whether we can say anything about $c$ _without_ solving the problem at all, the way the direct method bounded $\mathrm{max}_{(j)}|A^{n}_{j}|$ in {eq}`eq85` without ever integrating the scheme. We can, and the tool is the energy method.

The first step is to multiply {eq}`eq101` by $\tilde{v}^{*}$, the complex conjugate of $\tilde{v}$, and integrate across the domain. Using $\tilde{v}\rightarrow 0$ as $y\rightarrow\pm\infty$ and integrating the second-derivative term by parts gives

```{math}
:label: eq108
\int \left(|\tilde{v}_y|^{2}+k^{2}|\tilde{v}|^{2}\right)dy + \int \frac{\bar{U}_{yy}}{\bar{U}-c}|\tilde{v}|^{2}dy = 0
``` 

The first integral is real and strictly positive, so every bit of information about growth is hiding in the second one. Writing $\frac{1}{\bar{U}-c}=\frac{\bar{U}-c_r+ic_i}{|\bar{U}-c|^{2}}$ and taking the imaginary part of {eq}`eq108`,

```{math}
:label: eq109
c_i \int \frac{\bar{U}_{yy}}{|\bar{U}-c|^{2}}|\tilde{v}|^{2}dy = 0
``` 

Therefore, if the flow is unstable ($c_i \neq 0$), the integral must vanish on its own, which requires $\bar{U}_{yy}$ to **change sign** somewhere in the domain. This is _Rayleigh's inflection point theorem_: an inflection point in the background flow is a necessary (not sufficient) condition for instability. For our profile $\bar{U}_{yy}=\frac{\bar{U}}{L}[\delta(y+L)-\delta(y-L)]$, a pair of delta functions of opposite sign sitting exactly at the two edges of the shear zone, so the criterion is satisfied. Those two delta functions are precisely what the jump conditions {eq}`eq103` and {eq}`eq106` were encoding.

We can squeeze much more out of the same trick. Following Howard, substitute $\tilde{v}=(\bar{U}-c)G$ into {eq}`eq101`, which turns it into the neater form

```{math}
:label: eq110
\frac{d}{dy}\left[(\bar{U}-c)^{2}\frac{dG}{dy}\right]-k^{2}(\bar{U}-c)^{2}G = 0
``` 

Now multiply by $G^{*}$ and integrate by parts as before. All that survives is

```{math}
:label: eq111
\int (\bar{U}-c)^{2} Q \, dy = 0, \qquad Q \equiv |G_y|^{2}+k^{2}|G|^{2}\geq 0
``` 

The weight $Q$ is real and non-negative, which is what makes this so powerful. Since $(\bar{U}-c)^{2}=(\bar{U}-c_r)^{2}-c_i^{2}-2ic_i(\bar{U}-c_r)$, the imaginary part of {eq}`eq111` gives $c_i\int(\bar{U}-c_r)Q\,dy=0$, so for an unstable mode

```{math}
:label: eq112
c_r = \frac{\int \bar{U}Q \, dy}{\int Q \, dy}
``` 

i.e., the phase speed of a growing wave is a **weighted average of the background flow**, and so must lie between $U_{min}$ and $U_{max}$. This has an important physical consequence: there is always some $y_c$ with $\bar{U}(y_c)=c_r$. That location is called the _critical layer_, and it is where the wave exchanges momentum with the mean flow.

The real part of {eq}`eq111` gives $\int[(\bar{U}-c_r)^{2}-c_i^{2}]Q\,dy=0$. Combine it with the observation that $(\bar{U}-U_{min})(U_{max}-\bar{U})\geq 0$ everywhere, so that $\int(\bar{U}-U_{min})(U_{max}-\bar{U})Q\,dy\geq 0$, and after a few lines of algebra all the integrals cancel and we are left with _Howard's semicircle theorem_:

```{math}
:label: eq113
\left(c_r-\frac{U_{max}+U_{min}}{2}\right)^{2}+c_i^{2} \leq \left(\frac{U_{max}-U_{min}}{2}\right)^{2}
``` 

In words: for **any** parallel shear flow, the complex phase speed of any unstable mode must lie inside the semicircle in the upper half of the complex $c$ plane whose diameter is the range of the background velocity. We never solved an O.D.E. to get this. In particular the growth rate obeys $kc_i \leq k\frac{U_{max}-U_{min}}{2}$, which is the continuous analogue of the $|\lambda|\leq 1$ requirement from the von Neumann method: a bound on amplification obtained by inspection rather than by integration.

Let's check our own answer against it. For the profile in [FIG18](FIG18), $U_{max}=\bar{U}$ and $U_{min}=-\bar{U}$, so the semicircle is centred at $c_r=0$ with radius $\bar{U}$. In the unstable range the right hand side of {eq}`eq107` is negative, $c$ is purely imaginary, and $c_r=0$ sits exactly at the centre of the circle, which is forced by the antisymmetry of $\bar{U}$. The critical layer is then at $y=0$, the middle of the shear zone, which is just where we would expect the wave to feed on the mean flow. The amplitudes from {eq}`eq107` are:

| $\alpha L$ | $c_r/\bar{U}$ | $c_i/\bar{U}$ | inside semicircle? |
|---|---|---|---|
| $\rightarrow 0$ | 0 | $\rightarrow 1$ | touches the boundary |
| 0.1 | 0 | 0.871 | yes |
| 0.3984 | 0 | 0.505 | yes (fastest growing) |
| 0.5 | 0 | 0.368 | yes |
| 0.6392 | 0 | 0 | yes (neutral cutoff) |

Every unstable mode sits comfortably inside, and the bound is only saturated in the limit $\alpha L\rightarrow 0$. That limit is the vortex sheet, where the shear zone collapses to a discontinuity and we recover the classical Kelvin-Helmholtz result $c=\pm i\bar{U}$, sitting right at the top of the semicircle. Note that even though $c_i\rightarrow\bar{U}$ there, the growth rate $\sigma=kc_i$ still tends to zero because $k\rightarrow 0$; the actual maximum is $\sigma_{max}\approx 0.201\frac{\bar{U}}{L}$ at $\alpha L\approx 0.40$.

Finally, {eq}`eq107` tells us *why* there is a short-wave cutoff. The term $e^{-4\alpha L}$ is the square of the factor $e^{-2\alpha L}$ by which a disturbance at one edge ($y=+L$) is felt at the other ($y=-L$), a distance $2L$ away. Each edge on its own supports a neutral vorticity wave; instability is these two waves phase-locking and amplifying one another. Once $\alpha L$ is large the coupling $e^{-2\alpha L}$ becomes negligible, the two edge waves no longer feel each other, and the flow is stable. So Howard's theorem gives us the box that $c$ must live in, and {eq}`eq107` tells us whereabouts in that box this particular flow actually lands. 