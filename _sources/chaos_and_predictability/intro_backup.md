# Syllabus 

This graduate-level course will help you walk through some fundamental ideas in chaos (a.k.a. Butterfly effect), predictability and information theory. The goal is to learn how to (1) qualitatively and quantatively describe "predictability", (2) physically interpret a dynamical systems and (3) find the optimal patterns favoring model error growth in various systems. Prerequisites for this course including: applied math (ODE, PDE and linear algebra), numerical analysis (Python or Julia), statistics and atmospheric dynamics.   


:::{note}
This is the v1.0 handout of the course "Chaos and Predictability" taught in the [Department of Atmospheric Science](http://www.as.ntu.edu.tw/index.php/eng), [National Taiwan University](https://www.ntu.edu.tw/chinese2007/english/index.html). There might be tons of typos. If you find any error, feel free to DM me (email Dr. Kai-Chih Tseng: kaichiht@princeton.edu)
:::

## Course Outlines
__Part I: Lecturing__
* Week 1: __Predictability of Weather and Climate__
	* Climatological distribution vs forecast distribution 
    * State-dependent predictability and the mathematical assumptions 
    * Where the uncertainty comes from?
    * Ensemble forecast in weather and climate
    * HW1: Lorenz 63 model (testing state-dependent predictability) (10%)    
* Week 2: __Predictability Source in atmosphere__
	* The origins of predictability (Scale Seperation)
	* Brief introduction to balanced dynamics (Rossby wave dynamics)
	* Teleconnections 
* Week 3: __The minimalist models for studying chaos and predictability I: Lorenz 63__  
  	* The physical and mathematical background
    * Perfect model experitment
* Week 4: __The minimalist models for studying chaos and predictability II: Lorenz 96__
	* The physical and mathematical background
	* Single Scale
	* Multi-scale interations
	* Layapunov Exponent
	* The early and late stages of error growth in multi-scale flow. 
	* HW2: Lorenz 96 model (calculating the Layapynov exponent) (15%)
* Week 5: __The connection between statistical mechanics and ensemble forecasts I: theory__
	* Introduction to Liouville equation 
	* Liouville equation as a function of initial state vs current state 
	* Solution to Liouville equation
	* Connections between Liouville equation and ensemble forecasts
* Week 6: __The connection between statistical mechanics and ensemble forecasts II: applications__
	* Liouville equation in 1st-order ODE case
	* Liouville equation in Lorenz 63 model
	* Liouville equation and Singular Value Decomposition
	* Challengs of applying Liouville Equation to weather forecast
	* HW3: Derivation of Liouville equation + 1st-order ODE case simulation (15%)
* Week 7: __The stability theorem__
	* Stability in Geophysical Fluid Dynamics 
	* von Neumann method
	* Energy Method
* Week 8: __The generalized stability theorem__
	* Generalized stability theorem
	* Modal vs non-modal systems
    * Reduced-dimension problems
* Week 9: __Information theorem__
	* The connection between Layapunov exponent and perfect model experiments
	* Average Predictability Time
	* Predictable Components
	* HW4: Identifying the predictable components in a barotropic flow(20%)
* Week 10: __Ensemble forecast and data assimilation__
	* Ensemble forecast and data assimilation, the same problem or not?
	* Ensemble Kalman Filter 
* Week 11: __Future projection problems__
	* Initial value problems vs boundary value problems
	* Generalized stability in a non-linear/boundary value problem
	* Thompson et al. (2015) and Large-ensemble simulations 



__Part II: Literature Reviews__ (35%)
* Week 10: __The predictability of a flow which possesses many scales of motions (Lorenz 1969)__
	* Key Lecture: Week 8	
* Week 11: __Stochastic Forcing of ENSO by the Intraseasonal Oscillation (Moore and Kleeman 1999)__
	* Key Lecture: Week 9
* Week 12: __Ensemble-based sensitivity analysis (Hakim and Torn 2008)__
	* Key Lecture: Week 9
* Week 13: __The Critical Role of Non-Normality in Partitioning Tropical and Extratropical Contributions to PNA Growth (Henderson et al. 2020)__
	* Key Lecture: Week 9
* Week 14: __Physically Interpretable Neural Networks for the Geosciences: Applications to Earth System Variability (Toms et al. 2020)__
	* Key Lecture: Week 11

## Grading
HW X 4 (10%, 15%, 15% 20%)
Literature Review (35%)
Questions in class or office hour (5%)





	





 
And here is a code block:

```
e = mc^2
```

Check out the content pages bundled with this sample book to see more.
