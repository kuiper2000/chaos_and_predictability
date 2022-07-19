(week5)=
# Week 5: The connection between statistical mechanics and ensemble forecast I: Theory

In the following two weeks, we will go through some fundamental theory behind the *chaos and predictability*. As we know, the forecast uncertainty arises from two components (1) imperfect model physics and (2) imperfect observation in the initial condition. While most of the studies discussed from week 1 to week 4
were focusing on the intrinsic predictability limit (i.e., the upper limit of prediction), some efforts have been dedicated to the quantification/prediction of forecast error via statistical mechanics. In the context of predicting the forecast error, the Liouville equation arises as the general framework to describe the time evolution of ensemble density (i.e., the number of ensemble member in a finite phase space). 

## Introduction to Liouville equation
Liouville equation carries the name of a French mathematician, Joseph Liouville. You might have heard his name in the applied math course, when we introduced the Sturm-Liouville theory. and yes, these two "Liouville" are the same person. Liouville equation was first proposed by Joseph Liouville in 1838, which was about 100-year later after Euler proposed the Euler equation of fluid dynamics (1752,1757). While these two equations are highly similar, Liouville and Euler did carry out the formulas separately. Liouville equation (LE hereafter) can be considered as a more general form in Hamiltonian mechanics, which can be applied to not only fluid mechanics but also "any" dynamical system which is continuous in time and phase space. 

:::{note}
One should notice that we have a very strong assumption here, which is the continuity of a dynamical system. Whether such assumption holds in Navier-Stoke equation is still a million-dollar question (see [Millennium questions](https://www.claymath.org/millennium-problems) by Clay Institute of Mathematics . 
:::

## Jacobian Matrix and determinant
Before diving right in the Liouville equation, we will start with something simple, which will enable us to better understand the physical meaning behind. Two major ingredients in Liouville equation are (1) Jacobian matrix and (2) determinant. 

Mathematically, Jacobian matrix represents mapping one matrix to the other with a simple linear transformation. For example, when we manage to convert Polar coordinate to Cartesian coordinate, we can use the following formula, 

```{math}
:label: eq30
\begin{bmatrix}
x  \\
y   
\end{bmatrix} = 
\begin{bmatrix}
X(r,\theta) \\
Y(r,\theta) 
\end{bmatrix}
=
\begin{bmatrix}
r\mathrm{cos}(\theta)  \\
r\mathrm{sin}(\theta)   
\end{bmatrix}
```  
However, if the transformation is too complicated, we might want to find an easier way to approach the problem. One way to do that is through a simple linear transformation. So... how does it work? Let's start from a nonlinear function $y=f(x)$, which maps $x$ to $y$ or $f(x)$. If we only focus on a small range (let's say...centering around $p$) and assume a linear function $T(x)$ can well approximate the more complicated function $f(x)$, i.e., $T(x) = Ax+b$. This will give us 

```{math}
:label: eq31
T(x) = f(x) =Ax+b
```  

or 
```{math}
:label: eq32
T(x) = A(x-p)+f(p)
```  

One way to test if the linear assumption holds is by examining whether the $f(x)$ will approach $T(x)$ when $x\rightarrow p$, i.e.,

```{math}
:label: eq33
\mathbf{lim}_{x\rightarrow p} \frac{f(x)-T(p)-A(x-p)}{||x-p||} = 0
```  
{eq}`eq33` can be read as "if the difference between $f(x)$ and $T(x)$ is well approximated by linear function $A(x-p)$ and differentiable at point $p$".  

Now, we can extend the formula in {eq}`eq33` to higher dimensions. Let's first assume that $x$ has matrix form $\mathbf{x}$ with $n$ dimensions. In addition, we introduce another matrix $\mathbf{e}=[e_1,e_2,e_3...,e_j,...e_n]$, which is an n-dimension basis. We further assume that $\mathbf{e}=0$ except that its $j$th element equals 1 (i.e.,$e_{j}=1$) and include a small (nearly 0) constant, $h$. If $f(x)$ is differentiable at any given dimension, then we can rewrite {eq}`eq33` as 
```{math}
:label: eq34
\mathbf{lim}_{\mathbf{x}\rightarrow \mathbf{p}} \frac{f(\mathbf{p}+h\mathbf{e_j})-f(\mathbf{p})-A(h \mathbf{e}_j)}{||h\mathbf{e}_j||} = 0
```   

and because $A(h \mathbf{e_j})$ = $hA(\mathbf{e_j})$. {eq}`eq34` will lead to 
```{math}
:label: eq35
\mathbf{lim}_{h\rightarrow 0} \frac{f(\mathbf{p}+h\mathbf{e_j})-f(\mathbf{p})}{||h\mathbf{e}_j||} =A \mathbf{e}_j
```  

On the left-hand-side of {eq}`eq35` is the partial derivative of $f$ at point $\mathbf{p}$. Therefore, taking the derivative at $j$th dimension can be written as...     

```{math}
:label: eq36
\frac{\partial f}{\partial x_j} = 
\begin{bmatrix}
\frac{\partial f_1}{\partial x_j}  \\
\frac{\partial f_2}{\partial x_j}  \\
\frac{\partial f_3}{\partial x_j}  \\
\vdots \\
\frac{\partial f_m}{\partial x_j}   
\end{bmatrix}
```  

The full $\mathbf{A}$ matrix can be derived by taking the partial derivative of every component in $f$, i.e.,
```{math}
:label: eq37
\mathbf{A} = 
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots  & \frac{\partial f_1}{\partial x_j} \\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots  & \frac{\partial f_2}{\partial x_j}  \\
\frac{\partial f_3}{\partial x_1} & \frac{\partial f_3}{\partial x_2} & \cdots  & \frac{\partial f_3}{\partial x_j}  \\
\vdots \\
\frac{\partial f_m}{\partial x_1} & \frac{\partial f_m}{\partial x_2} & \cdots  & \frac{\partial f_m}{\partial x_j}
\end{bmatrix}
```  
 
Now, we can go back to our original question: converting the polar coordinate to Cartesian coordinate. We can find the $\mathbf{A}$ matrix for this problem is 

```{math}
:label: eq38
\mathbf{A} = 
\begin{bmatrix}
\frac{\partial r\mathrm{cos}(\theta)}{\partial r}  & \frac{\partial  r\mathrm{cos}(\theta)}{\partial \theta}\\
\frac{\partial r\mathrm{sin}(\theta)}{\partial r}  & \frac{\partial  r\mathrm{sin}(\theta)}{\partial \theta}
\end{bmatrix}
=
\begin{bmatrix}
\mathrm{cos(\theta)}   & -r\mathrm{sin(\theta)}\\
\mathrm{sin(\theta)}  & r\mathrm{cos(\theta)}
\end{bmatrix}
```  
Jacobian matrix is widely used in Geoscience and Geophysical Fluid Dynamics. One example is estimating the contribution of regional warming to the climate sensitivity, i.e., how much the earth should warm when we double the CO2 concentration. The role of Jacobian matrix in climate sensitivity studies is that we can consider the final coordinate as "the overall radiative forcing by greenhouse gases" and the initial coordinate (r.h.s) is the warming in SST over different regions. Thus, the yielded Jacobian matrix can tell us the relative importance of warming over different regions. In climate sensitivity study, it has an alternative name, Green's function. Some of the most updated studies about climate sensitivity can check this [workshop](https://usclivar.org/meetings/pattern-effect-workshop)

:::{note}
While pattern effect (or patching method) is a useful tool to identify the key regions responsible for global warming, its linear assumption is usually dubbed because the change in global temperature by CO2 forcing is not necessary small. Thus, there is a growing community turning their focus on the nonlinear effect of equilibrium climate sensitivity as well as it's local maximum (i.e., when global temperature is warm enough, the positive feedback from CO2 will get weaker (bell-like curve)) 
:::

The other key ingredient we need to learn is "determinant". Determinant evaluates how much the area spanned by the eigen basis has changed after the linear transformation. For example, in week 1, we have used the following case as an example, 

```{math}
:label: eq39
\begin{bmatrix}
1 & 0 \\
0 & 2 
\end{bmatrix} =
\begin{bmatrix}
1 & 0 \\
0 & 2 
\end{bmatrix} \begin{bmatrix}
1 & 0 \\
0 & 1 
\end{bmatrix}. 
```  
where the original basis, $[1,0]$ and $[0,1]$ (x-axis and y-axis respectively), is moved to $[1,0]$ and $[0,2]$ after transformation. We can easily find that the area originally spanned by $[1,0]$ and $[0,1]$ has grown twice bigger (from 1 $\rightarrow$ 4). The change in area size is equivalent to the determinant of the linear operator. i.e., 
$\mathrm{det}{(\begin{bmatrix}
1 & 0 \\
0 & 2 
\end{bmatrix})}$

Now, both ingredients are ready, the only question left is what is the connection among Jacobian matrix, determinant and Liouville equation? We can take a look of the example below. If today, we would like to convert a vector $\begin{bmatrix}
u  \\
v 
\end{bmatrix}$ to a vector $\begin{bmatrix}
x  \\
y 
\end{bmatrix}$
the corresponding linear transformation can be written as

```{math}
:label: eq40
\begin{bmatrix}
x  \\
y  
\end{bmatrix} =
\begin{bmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\
\frac{\partial y}{\partial u} & \frac{\partial x}{\partial v} 
\end{bmatrix} \begin{bmatrix}
u  \\
v  
\end{bmatrix}. 
```  

if we are using $du$ and $dv$ to represent the unit vectors in the original coordinate, then the unit vector in the final coordinates will be... 

```{math}
:label: eq41
\begin{bmatrix}
dx  \\
dy  
\end{bmatrix} =
\begin{bmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\
\frac{\partial y}{\partial u} & \frac{\partial x}{\partial v} 
\end{bmatrix} \begin{bmatrix}
du  \\
dv  
\end{bmatrix}. 
```  

for the area spanned by $dx$ and $dy$ in the final coordinate, it can be written as
```{math}
:label: eq42
\begin{align*}
& dxdy 
\\
& =||\mathrm{det}\begin{bmatrix}
dx & 0 \\
0  & dy  
\end{bmatrix}|| 
\\
&= ||\mathrm{det}
\begin{bmatrix}
\frac{\partial x}{\partial u}du & \frac{\partial x}{\partial v}dv \\
\frac{\partial y}{\partial u}du & \frac{\partial x}{\partial v}dv 
\end{bmatrix} || 
\\
&= 
||\mathrm{det}
\begin{bmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\
\frac{\partial y}{\partial u} & \frac{\partial x}{\partial v} 
\end{bmatrix} || dudv
\\
&=
||\mathrm{det} J(u,v)||dudv
\end{align*}

```   
from {eq}`eq42`, we can easily observe that the area spanned by unit vector has expanded (shrunk) by a factor of $||\mathrm{det}J(u,v)||$. 

## Connections between Liouville equation and ensemble forecasts

OK, now we know the mathematical meaning (geometric meaning) of a Jacobian matrix and determinant. To further understand the their roles in Liouville equation, let's first take a what Joseph Liouville found in 1838. For a dynamical system with third order differential terms  

```{math}
:label: eq43
x^{'''}=\mathbf{\Phi}(t,x,x^{''},x^{'''})
```  

he subsequently assumed that the solution has a form of  

```{math}
:label: eq44
x = x(t,a,b,c)
```  
where a, b and c are the initial state of this system. 
He further defined a Jacobian matrix, which maps the initial state to the final state, i.e., 

```{math}
:label: eq45
x(t=n_t,...) = J(a,b,c)x(t=0,a,b,c)
```  
the determinant of this Jacobian matrix can be written as

```{math}
:label: eq46
\rho(a,b,c,t)  
=||\mathrm{det}
\begin{bmatrix}
\frac{\partial x}{\partial a}      & \frac{\partial x}{\partial b} & \frac{\partial x}{\partial c} \\
\frac{\partial x^{'}}{\partial a}  & \frac{\partial x^{'}}{\partial b} & \frac{\partial x^{'}}{\partial c} \\  
\frac{\partial x^{''}}{\partial a} & \frac{\partial x^{''}}{\partial b} & \frac{\partial x^{''}}{\partial c}
\end{bmatrix} ||
```

He then proved that 
```{math}
:label: eq47
\frac{\partial \rho}{\partial t}=\rho\frac{\partial \mathbf{\Phi}}{\partial x^{''}}
```
:::{note}
{eq}`eq47` can be proved by using the chain rule i.e., 
$\frac{\partial x^{'''}}{\partial a}=\frac{\partial P}{\partial a} 
= \frac{\partial P}{\partial x}\frac{\partial x}{\partial a} + \frac{\partial P}{\partial x^{''}}\frac{\partial x^{''}}{\partial a}+\frac{\partial P}{\partial x^{'''}}\frac{\partial x^{'''}}{\partial a}$  (HW)
:::

Let's take a look what {eq}`eq47` tells us. First, {eq}`eq44` describes how the dynamical system gonna evolve with time. If $x$ represents the forecast uncertainty or forecast error (i.e., the distance between any given ensemble member to the ensemble mean), then {eq}`eq44` is a "prognostic equation of model uncertainty". Second, the determinant tells us the how the volume spanned by eigen basis changes after the linear transformation. Physically, the volume spanned by different ensemble members is so-called ensemble spread. Thus, {eq}`eq47` is a "prognostic equation for ensemble spread" (i.e., forecast the entire PDF). In addition, it also tells us the time tendency of ensemble spread is only related to the highest-order differentiation i.e., $\frac{\partial }{\partial x^{'''}}$.  

{eq}`eq47` also has some very interesting applications in the information theory. If we rewrite {eq}`eq47` by dividing both side with $\rho$, we can have
```{math}
:label: eq48
\frac{\partial \mathrm{ln}\rho}{\partial t}=\frac{\partial \mathbf{\Phi}}{\partial x^{''}}
```
we further integrate both sides over the entire phase space. One can find that the integration of right-hand-side of {eq}`eq49` is always 0 (why?). Thus, {eq}`eq48` is also equivalent to 

```{math}
:label: eq49
\int\frac{\partial \mathrm{ln}\rho}{\partial t}dx^{''}=\int\frac{\partial \mathbf{\Phi}}{\partial x^{''}} dx^{''}=0
```
this also gives us $\frac{\partial }{\partial t} \int\mathrm{ln}\rho dx^{''}$ = 0 (note that we swap the order of integration and $\frac{\partial }{\partial t}$ since both of them are linear operators). {eq}`eq49` is so-called the "conservation of information". In [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)#Definition), the definition of "information" is very similar to the left-hand-side of {eq}`eq49`, where the information provided by $i$th element is $\mathrm{ln}\rho(x_{i})$.


Now, we can imagine how powerful {eq}`eq47` is since it forecasts the entire probability density function. Different from the traditional probabilistic forecast method, which estimates the ensemble spread by generating tens of hundreds of ensemble simulations, the Liouville equation is a more powerful tool. However, there is no free meal. Implementing Liouville equation in real operational forecast is extremely challenging since we need to find the analytical solution of a dynamical system first (i.e., {eq}`eq44`) and then calculate its derivative over all possible dimensions (i.e., the right-hand-side of {eq}`eq47`). This is only possible when the dynamical system is simple enough. 


Next week, we will walk through a few cases where we can analytical implement the Liouville equation in the forecast of the PDF. 
:::{note}
The Nobel prize winner, Roger Penrose (Penrose (1989)), who was also Stephan Hawkins' Ph.D., committee, referred to Liouville equation as a very beautiful formula since it describes that the volume of any region of phase spaces must remain the same. Penrose also proposed that a black hole can destroy information (i.e., $\rho$ in {eq}`eq47` is not conserved) and thus the statistical mechanics might break down at the space-time singularity. 


In addition to his achievement in theoretical physics, Penrose is also an artist who is famous for his 3D Penrose tiling. An award-winning app (mobile game), [monument valley](https://en.wikipedia.org/wiki/Monument_Valley_(video_game)) is inspired by his arts. NYC has a National Museum of Mathematics where you can find some of Penrose's work.  
:::

```{figure} https://upload.wikimedia.org/wikipedia/commons/thumb/5/52/RogerPenroseTileTAMU2010.jpg/220px-RogerPenroseTileTAMU2010.jpg
---
name: FIG12
---
Roger Penrose and Penrose tiling
```


## References
```{bibliography} ../references.bib
:filter: docname in docnames
```

