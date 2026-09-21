(week3)=
# Week 3: The minimalist models for studying chaos and predictability I: Lorenz 63

In weeks 1 and 2, we learned how to define predictability and explore the potential predictability sources in Earth atmosphere. This week, we will talk about the origins of the studies of predictability and chaos. 

Next to relativity and statistical mechanics, chaos is probably one of the most groundbreaking discoveries in theoretical physics in 20th century. Although Edward Lorenz is not "the" father of chaos theory (some fundamental work had been done by [Henri Poincare](https://en.wikipedia.org/wiki/Henri_Poincar%C3%A9)), he was definitely one of a few pioneers who systematically unraveled this field. This week, we will walk through the details of his legendary 63 paper: "Deterministic non-periodic flow"{cite}`lorenz1963deterministic`. 


## The origins of the story

Some of you might have heard of the story that Edward Lorenz (Ed hereafter) accidentally found "chaos" when he was dealing with the rounding error problem. Specifically, when he tried to save the simulation output for the next experiments (while he was taking a break for a cup of coffee), he unpurposefully truncated the data which led a totally different simulated result at the end. He then realized that for a non-linear chaotic system, even an infinitesimal error in the initial state can ultimately lead to an unpredictable future i.e., the relation between aperiodicity and the loss of prediction skills.   

To better understand how the small initial error led to an chaotic future, Ed needed a system, which was both chaotic but also simple enough for understanding (also see this [workshop](https://eapsweb.mit.edu/news/2018/celebration-two-pioneers-modern-meteorology)  in honor of Lorenz and Charney). At the very beginning, Ed was working on a set of equation with 12 variables mimicking the state-of-art forecast system 

```{math}
:label: eq16
\begin{align*}
\frac{dX_i}{dt}=\sum_{m,n}a_{imn}X_mX_n+\sum_{m}b_{im}X_m+c
\end{align*}
```  
With the chosen parameters, he found some non-periodic behaviors for this set of equations. However, it was too complicated to boil down (and later in his career, he found it's impossible to simplify this equation). So...Ed continued searching for a much simpler system. His effort bore fruits when he visited Barry Saltzman in 1961 {cite}`saltzman1962finite`. Barry showed him 7 equations representing the convective system, where most of them showed periodic behavior except one refused to settle down. This was the famous Rayleigh-Bernard convection{cite}`rayleigh1916lix`, also the starting point of the famous Lorenz 63 model. 

:::{note}
Equation {eq}`eq16` is so-called Lorenz 96 model{cite}`lorenz1996predictability`, which describes non-linear waves propagating in a single zonal band. More details will be provided next week. One reason making it impossible to simplify is that when we try to solve the linearnized {eq}`eq16` (stability test, see week 7), it will lead a polynomial equation with orders higher than 5 (i.e., $\sum_{n=1}^{N} a_n x^n$ where N$\geq 5$. According to [Galois theory](https://en.wikipedia.org/wiki/Galois_theory), there is no analytical solution for a polynomial equation with order higher than 5 (or can't be represented with the existing math symbols such as +,-, X, / or $()^{\frac{1}{n}}$ ($\forall n\in N$) without introducing the concept of infinity).  
:::

## Rayleigh-Barnard Convection 

Rayleigh-Barnard problem is a set of equation describing the convective-driven/diffusion-driven flow. Obviously, there are two competing forces of mass transport in this kind of flow (1) buoyancy and (2) drag. So we need momentum and thermodynamics equations in Rayleigh-Barnard problem. 
```{math}
:label: eq17
\begin{align*}
\frac{\partial \mathbf{u}}{\partial t}+(\mathbf{u}\cdot\nabla)\mathbf{u} &=\nu\Delta\mathbf{u}-\frac{\nabla_h p}{\rho_0}+(\alpha T)g \\
\frac{\partial T}{\partial t}-u_z\Gamma+ (\mathbf{u}\cdot\nabla_h)T &=\kappa\Delta T \\
\nabla\cdot \mathbf{u} &= 0
\end{align*}
```
where $\mathbf{u}$ is 3-dimensional momentum in vector form, g is the gravitational acceleration, T is the departure of temperature from a reference temperature profile (or potential temperature since we will introduce Boussinesq approximation later). $\alpha$ is the coefficient of thermal expansion (not specific volume), which is define as $\frac{1}{A}\frac{dA}{dT}$ (A is volume and T is the temperature). The last term in first equation of {eq}`eq17` is so-called reduced gravity or buoyancy. $u_z$ is the vertical velocity and $\Gamma$ is the lapse rate, defined as $-\frac{(T_L-T_H)}{d}$ (d is the vertical depth of fluid) and $\kappa\Delta T$ is the thermal diffusion. The last equation of {eq}`eq17` is the mass continuity with Boussinesq approximation. To simplify the problem, we further implement two boundary conditions: (1) $\mathbf{u}$ = 0 at z=0,d and (2) $T$ = 0 at z=0,d. In addition, we assume $\mathbf{u}=(u_x,0,u_z$). These implemented conditions simplify the Rayleigh-Benard convection to problem of a bounded (at d=0 and z), 2D convection. 

We further introduce a few non-dimensional parameters to rescale the equation (,which can be done with Buckingham $\pi$ theorem as well).  
```{math}
:label: eq18
\begin{align*}
x^* &= \frac{x}{d} \\
z^* &= \frac{z}{d} \\
t^* &= \frac{t}{(d^2/\kappa)} \\
u_x^* &= \frac{u_x}{\kappa/d} \\
u_z^* &= \frac{u_z}{\kappa/d} \\
T^*   &= \frac{T}{(\kappa\nu/g\alpha d^3)} \\
\end{align*}
```
where $\kappa/d$ is the characteristic velocity and $t_c$ is the characteristic timescales defined as $(d^2/\kappa)$ (dissipation timescale). With this non-dimensional scaling, {eq}`eq17` is reduced to a set of equations with only two free parameters. In addition, since the solution is 2D and bounded, we can use mass continuity to find the relation between mass stream function and momentum, i.e., $\mathbf{u}=(u_x,0,u_z)=(-\frac{\partial \psi}{\partial z},0,\frac{\partial \psi}{\partial x})$.  This one-on-one relation is so-called invertibility, which is similar to the concept of PV inversion (i.e., {eq}`eq11`).  

The x component of {eq}`eq17` contains $\partial_x p$, and the z component contains $\partial_z p$. Hence by differentiating the former with respect to z and the latter with respect to x, we have two equations containing a $\partial_x\partial_z$ term, which drops out on subtraction. Substituting the stream function into this equation results in   

```{math}
:label: eq19
\begin{align*}
\partial_{t^*} \Delta^* \psi^*-\partial_{z^*} \psi^* \partial_{x^*} \Delta^* \psi^*+\partial_{x^*} \psi^*\partial_{z^*}\Delta^* \psi^* &= \sigma\Delta^{*2} \psi^*+\sigma\partial_{x^*} T^* \\
\partial_{t^*} T^* -\partial_{z^*} \psi^* \partial_{x^*} T^*+\partial_{x^*} \psi^*\partial_{z^*}T^* &= \Delta^* T^* +R \partial_{x^*} \psi^* 
\end{align*}
```

where $\sigma$ is Prandtl number, defined as $\frac{\nu}{\kappa}$, and $R$ is Rayleigh number, defined as $\frac{g\alpha d^3 (T_H-T_L)}{\kappa\nu}$. By varying $\sigma$ and $R$, we can acquire the statistics of Rayleigh-Bernard convection in different parameter spaces. The animation below is based on a simulation of $(\sigma/R)=\frac{1}{10^8}$ . 

```{figure} RB2D/animation/image.gif
---
name: FIG6
---
An example of buoyancy-driven Rayleigh-Barnard convection with $(\sigma/R)=\frac{1}{10^8}$ and an aspect ration of 16:9. (favoring the generation of zonal wave number 2 convection). The shading shows the temperature anomaly.  
```

````{admonition} Where does the Rayleigh number come from? (Reynolds number + terminal velocity)
:class: dropdown, tip

**Why the Reynolds number is not enough.** The Reynolds number is the ratio of inertial to viscous force on a parcel of size $L$ moving at speed $U$,

```{math}
\mathrm{Re}=\frac{\mathrm{inertia}}{\mathrm{viscous}}\sim\frac{\rho_0 U^2/L}{\mu U/L^2}=\frac{UL}{\nu}
```

The catch is that $U$ must be *supplied from outside* (a pipe, a fan, a jet). In Rayleigh-Barnard convection there is no imposed velocity at all: the only inputs are $g$, $\alpha$, $\Delta T=T_H-T_L$, $d$, $\nu$ and $\kappa$. So the strategy is to **build a velocity scale out of the buoyancy forcing itself, and then feed it into a Reynolds-type number**.

**The terminal velocity of a buoyant blob.** Take a warm parcel of size $d$ with temperature excess $\Delta T$. With $\alpha=\frac{1}{A}\frac{dA}{dT}$ its density deficit is $\delta\rho/\rho_0=-\alpha\Delta T$, so the buoyancy force per unit volume is $\rho_0 g\alpha\Delta T$ -- exactly the reduced-gravity term $(\alpha T)g$ in {eq}`eq17`. The parcel does *not* accelerate forever, because viscous drag resists it. For a slow (low-$\mathrm{Re}$) blob the drag per unit volume is $\sim\mu U/d^2=\rho_0\nu U/d^2$. Balancing buoyancy against drag gives the **terminal (Stokes) velocity**

```{math}
:label: eq19a
U_t \sim \frac{g\alpha\Delta T\, d^2}{\nu}
```

(exactly $U_t=\frac{2}{9}g\alpha\Delta T a^2/\nu$ for a sphere of radius $a$, from $F_b=\frac{4}{3}\pi a^3\rho_0 g\alpha\Delta T$ against Stokes drag $F_d=6\pi\mu a U$). This is the intrinsic velocity scale of the problem: nothing external was needed.

**First attempt.** Plugging {eq}`eq19a` into the Reynolds number gives

```{math}
\mathrm{Re}=\frac{U_t d}{\nu}=\frac{g\alpha\Delta T\, d^3}{\nu^2}\equiv \mathrm{Gr}
```

the Grashof number. This is *almost* the Rayleigh number, but it only counts one of the two ways a buoyant parcel can be defeated: momentum diffusion.

**The missing ingredient.** A rising warm blob has a second way to die -- it can lose its buoyancy before it gets anywhere, by conducting heat to its surroundings at rate $\kappa$. The blob only matters if it crosses the layer *before* its thermal anomaly is smeared out. So the resistance that belongs in the denominator is the **thermal** diffusivity, i.e. we want a Peclet number rather than a Reynolds number:

```{math}
\mathrm{Pe}=\frac{U_t d}{\kappa}=\frac{1}{\kappa}\cdot\frac{g\alpha\Delta T d^2}{\nu}\cdot d=\frac{g\alpha\Delta T\, d^3}{\kappa\nu}\;\equiv\;R
```

which is precisely the Rayleigh number. Three equivalent readings:

1. **As $\mathrm{Re}\times\mathrm{Pr}$.** Since $\sigma=\nu/\kappa$, $R=\mathrm{Gr}\cdot\sigma$. The Rayleigh number *is* a Reynolds number built on the terminal velocity, multiplied by the Prandtl number -- the two free parameters of {eq}`eq19` are just $\{\mathrm{Re},\mathrm{Pr}\}$ in disguise.
2. **As a ratio of timescales.** With $\tau_\kappa=d^2/\kappa$ (time to diffuse the anomaly away) and $\tau_{\mathrm{rise}}=d/U_t=\nu/(g\alpha\Delta T d)$ (time to cross the layer), $R=\tau_\kappa/\tau_{\mathrm{rise}}$. If $R\gg1$ the parcel crosses before it can be diffused away and we get convection; if $R\ll1$ it never makes it and we get pure conduction. This is also why $t_c=d^2/\kappa$ was the natural choice in {eq}`eq18`.
3. **As stratification times two dissipation times.** With $N^2=g\alpha\Delta T/d$, $R=N^2\tau_\nu\tau_\kappa$ where $\tau_\nu=d^2/\nu$: (driving strength) $\times$ (how long momentum survives) $\times$ (how long heat survives).
````

## Lorenz 63 model I (solutions) 

To see if equation {eq}`eq19` can be boiled down to the simplest form but still maintains its chaotic nature, we apply Fourier decomposition to {eq}`eq19` to simplify the problem. Since the solution is periodic in zonal direction but bounded in vertical direction, we assume it has the following form. 


```{math}
:label: eq20
\begin{align*}
T^*(x^*,z^*)    & = \sum_{k1=-\infty}^{\infty}\sum_{k2=1}^{\infty} \hat{T}_{k1,k2}e^{ik_1\pi x/l} \mathrm{sin}(k_2\pi z/d)\\
\psi^*(x^*,z^*) & = \sum_{k1=-\infty}^{\infty}\sum_{k2=1}^{\infty} \hat{\psi}_{k1,k2}e^{ik_1\pi x/l} \mathrm{sin}(k_2\pi z/d)
\end{align*}
```

Before we dive right in {eq}`eq19`, a simple math trick we should know is that the interaction of flow with different scales (i.e., wave number) can yield new scales or new wave number. i.e., $\psi_{k_1} \psi_{k_2} = \hat{\psi_{k_1}}\hat{\psi_{k_2}}e^{k_1 x}e^{k_2 x}=\hat{\psi_{k_1}}\hat{\psi_{k_1}} e^{(k_1+k_2) x}$. Therefore, to predict $\psi_{k_1}$ or $T_{k_1}$ in of {eq}`eq19`, we need at least one source term (e.g., advection) which has the same wave number after the nonlinear interaction.  Otherwise, we will never be able to predict $\psi_{k_1}$ or $T_{k_1}$ due to the orthogonality of Fourier modes (this also leads to another very interesting math problem: [Partition](https://en.wikipedia.org/wiki/Partition_(number_theory)), where [Srinivasa Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan) made some groundbreaking contributions.)  

Back to {eq}`eq19`, by inspecting {eq}`eq19`, the minimum way of sustaining the nonlinearity is to keep the advection term with the least amount of wave numbers. Here we choose the imaginary part of $\psi_{k_1=1}$ as the starting point. If we drop the advection term in the prognostic equation of $\psi$ to simplify the equation, the only way to keep its nonlinearity is through $T$ (since $\Delta^2 \psi$ is a damping term). Thus, we must have a prognostic equation for $\hat{T_{k1=1}}$. On the other hand, due to the form of $\sigma\partial_{x}T$ in {eq}`eq19`, the prognostic equation for $\hat{T_{k1=1}}$ must predict its real part (I will let you think about why). 

Because we have dropped the advection term in the prognostic equation of $\psi$, we must keep the advection term in the prognostic equation of $T$ to ensure the nonlinearity of the whole dynamical system. According to our previous discussion, as long as the difference in wave number between two scales equals 1 i.e., $k_1+k_2=1$ or $k_1-k_2=1$, the nonlinear interaction between $k_1$ and $k_2$ can always generate wave number 1 flow. Thus, the simplest way to include the forcing term for $\hat{T_{k_1=1}}$ without introducing too many new scales is to have $\hat{T_{k_1=0}}$ since we already have $\hat{\psi_{k_1=1}}$. (i.e., wave number 1 wind field interacts with the mean temperature). This will bring us the third variable $\hat{T_{k_1=0}}$. 

Let's take a quick look at what we have so far: the real part of $\hat{\psi_{k1=1}}$, the imaginary part of $\hat{T_{k_1=1}}$ and $\hat{T_{k_1=0}}$. To see if we can close the loop (i.e., enough forcing terms for all included wave numbers), we can substitute $T$ in {eq}`eq19` with $\hat{T_{k_1=0}}$ and check if there is any term can be used as the forcing to generate $\hat{T_{k_1=0}}$. It turned the nonlinear interaction between $\hat{\psi_{k_1=1}}$ and $\hat{T_{k_1=1}}$ can yield $\hat{T_{k_1=0}}$. Thus, we have officially closed the loop!! (Woohoo)

Repeat the same practice, we can find the necessary meridional wave numbers (I will let readers explore this part). Therefore, the solution of {eq}`eq19` is boiled down to 

```{math}
:label: eq21
\begin{align*}
\psi^*(x^*,z^*)    & = X(t)\mathrm{sin}(\pi ax^*)\mathrm{sin}(\pi z^*) \\
T^*(x^*,z^*)       & = Y(t)\mathrm{cos}(\pi ax^*)\mathrm{sin}(\pi z^*)-Z(t)\mathrm{sin}(2\pi z^*)
\end{align*}
```  
where $a$ is the aspect ratio of domain (i.e., d/l) and $X(t)$, $Y(t)$ and $Z(t)$ correspond to the real part of $\hat{\psi_{k_1=1,k_2=1}}$, the imaginary part of $\hat{T_{k_1=1,k_2=1}}$ and $\hat{T_{k_1=0,k_2=2}}$ respectively. By substituting $T$ and $\psi$ in {eq}`eq19` with the definition in {eq}`eq21` and removing the constant. We can have the famous Lorenz 63 model 

```{math}
:label: eq22
\begin{align*}
\dot{X} &= -\sigma X+ \sigma Y \\
\dot{Y} &= -XZ+ rX- Y \\
\dot{Z} &= XY- bZ \\
\end{align*}
```  
where $\sigma$ is Prandtl number, defined as $\nu/\kappa$, r is the ratio between $R$ (Rayleigh number) and $Rc$ (critical Rayleigh number, defined as $\pi^4 a^{-2}(1+a^2)^{3}$) and $b$ is $4(1+a^2)^{-1}$ (see more details in {cite}`rayleigh1916lix` and {cite}`saltzman1962finite`). $Rc$ provides a mean of estimating if the flow is dominated by buoyancy forcing or turbulent diffusion. When $R$ is much bigger than $Rc$, the diffusion is small compared to other terms in {eq}`eq16`, thus the small-scale features will take longer to decay (this is the reason why we are more likely to see a turbulent system). In a domain with $a^2=\frac{1}{2}$, the critical value of $R$ happens at $\frac{27\pi^4}{4}$. 

````{admonition} The full derivation, from the Boussinesq equations (17) to Lorenz 63 (22)
:class: dropdown, tip

There are four steps: **(A)** eliminate pressure to get the vorticity-streamfunction form; **(B)** non-dimensionalise to get {eq}`eq19`; **(C)** Galerkin-truncate with {eq}`eq21`; **(D)** rescale amplitudes and time to get {eq}`eq22`.

**Step A -- eliminate the pressure ({eq}`eq17` $\rightarrow$ vorticity form).** With $\mathbf{u}=(-\partial_z\psi,0,\partial_x\psi)$, continuity is satisfied identically, $\partial_x u_x+\partial_z u_z=-\partial_x\partial_z\psi+\partial_z\partial_x\psi=0$. Writing the two momentum components,

```{math}
\begin{align*}
\partial_t u_x+(\mathbf{u}\cdot\nabla)u_x &= \nu\Delta u_x-\tfrac{1}{\rho_0}\partial_x p \\
\partial_t u_z+(\mathbf{u}\cdot\nabla)u_z &= \nu\Delta u_z-\tfrac{1}{\rho_0}\partial_z p+g\alpha T
\end{align*}
```

apply $\partial_z$ to the first, $\partial_x$ to the second, and subtract. Both pressure terms become $\partial_x\partial_z p$ and cancel. Defining the only surviving vorticity component

```{math}
\zeta\equiv\partial_x u_z-\partial_z u_x=\partial_{xx}\psi+\partial_{zz}\psi=\Delta\psi
```

and noting $\mathbf{u}\cdot\nabla=-\partial_z\psi\,\partial_x+\partial_x\psi\,\partial_z$ (the Jacobian $J(\psi,\cdot)$), we get

```{math}
\begin{align*}
\partial_t\Delta\psi-\partial_z\psi\,\partial_x\Delta\psi+\partial_x\psi\,\partial_z\Delta\psi &= \nu\Delta^2\psi+g\alpha\,\partial_x T \\
\partial_t T-\partial_z\psi\,\partial_x T+\partial_x\psi\,\partial_z T &= \kappa\Delta T+\Gamma\,\partial_x\psi
\end{align*}
```

where the last term used $u_z=\partial_x\psi$ and $\Gamma=(T_H-T_L)/d$.

**Step B -- non-dimensionalise ($\rightarrow$ {eq}`eq19`).** Using {eq}`eq18` plus $\psi=\kappa\psi^*$ (since $[\psi]=$ velocity $\times$ length $=\kappa$): in the vorticity equation every term on the left scales as $\kappa^2/d^4$, while the right-hand terms are $\nu\kappa d^{-4}\Delta^{*2}\psi^*$ and $g\alpha\cdot d^{-1}\cdot\frac{\kappa\nu}{g\alpha d^3}\partial_{x^*}T^*=\kappa\nu d^{-4}\partial_{x^*}T^*$. Dividing through by $\kappa^2/d^4$ leaves the single parameter $\sigma=\nu/\kappa$ in front of *both* right-hand terms -- which is exactly why $T_c=\kappa\nu/(g\alpha d^3)$ was chosen in {eq}`eq18`.

In the temperature equation all terms scale as $\kappa^2\nu/(g\alpha d^5)$ *except* the lapse-rate term, whose ratio to that scale is

```{math}
\Gamma\frac{\kappa}{d}\Big/\frac{\kappa^2\nu}{g\alpha d^5}=\frac{g\alpha d^4\Gamma}{\kappa\nu}=\frac{g\alpha d^3\Delta T}{\kappa\nu}=R
```

This is where the Rayleigh number is *born*: it is not inserted by hand, it is the only free parameter the heat equation can support. The result is {eq}`eq19`.

**Step C -- Galerkin truncation.** Drop the stars and write $k\equiv\pi a$, $m\equiv\pi$, $K^2\equiv k^2+m^2=\pi^2(1+a^2)$, so that {eq}`eq21` reads $\psi=X\sin kx\sin mz$ and $T=Y\cos kx\sin mz-Z\sin 2mz$.

For the vorticity equation, a single mode satisfies $\Delta\psi=-K^2\psi$, so the Jacobian vanishes **identically**:

```{math}
-\partial_z\psi\,\partial_x\Delta\psi+\partial_x\psi\,\partial_z\Delta\psi=-K^2\left(-\partial_z\psi\,\partial_x\psi+\partial_x\psi\,\partial_z\psi\right)=0
```

This is the rigorous version of "we drop the advection term in the $\psi$ equation" above -- it is not an approximation, it is exact for one mode. What remains is $-K^2\dot{X}=\sigma K^4X-\sigma kY$, i.e.

```{math}
:label: eq21a
\dot{X}=-\sigma K^2X+\frac{\sigma k}{K^2}Y
```

For the temperature equation the linear pieces are $\partial_t T=\dot{Y}\cos kx\sin mz-\dot{Z}\sin 2mz$, $\Delta T=-K^2Y\cos kx\sin mz+4m^2Z\sin 2mz$ and $R\partial_x\psi=RkX\cos kx\sin mz$. The nonlinear pieces, using $\psi_z=Xm\sin kx\cos mz$ and $\psi_x=Xk\cos kx\sin mz$, are

```{math}
\begin{align*}
-\psi_z T_x &= kmXY\sin^2\!kx\,\sin mz\cos mz \\
\psi_x T_z  &= kmXY\cos^2\!kx\,\sin mz\cos mz-2kmXZ\cos kx\,\sin mz\cos 2mz
\end{align*}
```

Adding them, $\sin^2+\cos^2=1$ and $\sin mz\cos mz=\frac{1}{2}\sin 2mz$, so the first pair collapses to $\frac{km}{2}XY\sin 2mz$ -- the source that feeds $Z$. For the last term, $\sin mz\cos 2mz=\frac{1}{2}[\sin 3mz-\sin mz]$, **and here is the crucial truncation**: $\sin 3mz$ lies outside the retained basis and is discarded (this is the "nonlinear interaction generates new wave numbers" argument above, made concrete). The surviving $-\frac{1}{2}\sin mz$ returns $+kmXZ\cos kx\sin mz$. Projecting onto $\cos kx\sin mz$ and onto $\sin 2mz$:

```{math}
:label: eq21b
\begin{align*}
\dot{Y} &= RkX-K^2Y-km\,XZ \\
\dot{Z} &= \tfrac{km}{2}XY-4m^2Z
\end{align*}
```

The Lorenz structure is already visible in {eq}`eq21a`-{eq}`eq21b`: linear coupling $X\leftrightarrow Y$, the $-XZ$ feedback, the $+XY$ source for $Z$, and three distinct damping rates.

**Step D -- rescale.** Let $\tau=K^2t$ and

```{math}
X=\frac{\sqrt{2}K^2}{km}\hat{X},\qquad Y=\frac{\sqrt{2}K^6}{k^2m}\hat{Y},\qquad Z=\frac{K^6}{k^2m}\hat{Z}
```

Substituting and dropping hats gives exactly {eq}`eq22`, with the parameters falling out automatically:

```{math}
r=\frac{Rk^2}{K^6}=\frac{R\pi^2a^2}{\pi^6(1+a^2)^3}=\frac{R}{R_c},\qquad R_c=\pi^4a^{-2}(1+a^2)^3,\qquad b=\frac{4m^2}{K^2}=\frac{4}{1+a^2}
```

Equivalently, one can state {eq}`eq21` in already-normalised form (Saltzman's own convention) and go straight to {eq}`eq22` with no rescaling step:

```{math}
\begin{align*}
\psi^*(x^*,z^*) &= \frac{\sqrt{2}(1+a^2)}{a}X(\tau)\,\mathrm{sin}(\pi ax^*)\mathrm{sin}(\pi z^*) \\
T^*(x^*,z^*)    &= \frac{R_c}{\pi}\left[\sqrt{2}\,Y(\tau)\,\mathrm{cos}(\pi ax^*)\mathrm{sin}(\pi z^*)-Z(\tau)\,\mathrm{sin}(2\pi z^*)\right]
\end{align*}
```

with $\tau=\pi^2(1+a^2)t^*$.

**Physical dictionary.** $X$ is $\mathrm{Re}\,\hat{\psi}_{1,1}$, the intensity of the convective roll; $Y$ is $\mathrm{Im}\,\hat{T}_{1,1}$, the horizontal temperature contrast between the rising and sinking branches; $Z$ is $\hat{T}_{0,2}$, the departure of the mean vertical profile from pure conduction. The term $-XZ$ in $\dot{Y}$ says that as the roll flattens the mean profile it starves its own driving -- this is the nonlinearity responsible for chaos. The term $+XY$ in $\dot{Z}$ says the roll transports heat and flattens the profile. The three damping rates $\sigma$, $1$ and $b$ are simply $\nu K^2$, $\kappa K^2$ and $\kappa(2m)^2$ after the rescaling.
````

From {eq}`eq19` to {eq}`eq22`, we can find the Lorenz 63 originated from a heavily truncated Rayleigh-Barnard problem, where the three variables, $X(t),Y(t)$ and $Z(t)$ represent the first three Fourier modes.  

## Lorenz 63 model II (stability analysis)

An easy way to visualize how the chosen parameters affect the dynamical behavior in Lorenz 63 model is through the stability analysis. While some basic ideas are provided here, we will cover more details in week 7. The basic concept of linear stability analysis is examining the linear error growth rate and testing if the system is linearly stable or unstable (e.g., $\psi \sim e^{\lambda t}$, where $\lambda>0$ for unstable case and vice versa). If the given system is linearly unstable, it will deviate from linear regimes and become unpredictable in finite simulation time.

The linearized Lorenz 63 model can be written as:

```{math}
:label: eq23
\begin{bmatrix}
\dot{X}  \\
\dot{Y}  \\
\dot{Z} 
\end{bmatrix} = 
\begin{bmatrix}
-\sigma & \sigma & 0 \\
r & -1 & 0 \\
0 & 0 & -b
\end{bmatrix}
\begin{bmatrix}
X  \\
Y  \\
Z 
\end{bmatrix}
=\mathbf{L}\begin{bmatrix}
X  \\
Y  \\
Z 
\end{bmatrix}
```  
It's not hard to find that $Z(t)$ is now decoupled with $X(t)$ and $Y(t)$. Interestingly, $Z(t)$ also represents mean state (wave number 0) in the original Rayleigh-Barnard equation. Physically, linearization also assumes that the mean state doesn't interact with the transient states. 

By calculating the eigen value of $\mathbf{L}$ in {eq}`eq23`, we can estimate the growth rate/frequency of linearized Lorenz 63 model. Since we are looking for a chaotic system, that means we are interested in the case when $\lambda>0$ (roots are not only real but positive) in {eq}`eq24`. 

```{math}
:label: eq24
[\lambda+b][\lambda^2+(\sigma+1)\lambda+\sigma(1-r)]=0
```  
{eq}`eq24` has three real roots when $r>0$. One of them is always $\lambda=-b<0$; the other two are the roots of $\lambda^2+(\sigma+1)\lambda+\sigma(1-r)$, whose product is $\sigma(1-r)$. When $r<1$ that product is positive and both roots are negative (the conduction state is stable), whereas when $r>1$ the product turns negative, so there is exactly **one positive root** and the conduction state loses stability along a single unstable direction. This is consistent with the finding in Rayleigh 1912 {cite}`rayleigh1916lix` (i.e., the convection only happens when $R>Rc$). Another way to visualize the linear instability on a phase space is by estimating the divergence of physical states or the trace of the linear operator in {eq}`eq23`. i.e., $\frac{\partial\dot{X}}{\partial X}+\frac{\partial\dot{Y}}{\partial Y}+\frac{\partial\dot{Z}}{\partial Z}=\mathrm{Tr}(\mathbf{L})=-(\sigma+b+1)$ (see discussion in {ref}`week1`). If $\mathrm{Tr}(\mathbf{L})>0$, it indicates the volume spanned by different ensemble members will increase with the increase of forecast lead time. However, a value of $\mathrm{Tr}(\mathbf{L})<0$ doesn't mean this system will shrink to a single point. Instead, it may simply flatten to a plate. In most cases, {eq}`eq24` is not valid due to the non-stationary basic states. A more reasonable way to linearize the Lorenz 63 model will be discussed in week 5 when we talk about its connections with statistical mechanics.  


## References
```{bibliography} ../references.bib
:filter: docname in docnames
```





