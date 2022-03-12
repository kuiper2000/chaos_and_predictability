# Syllabus 

This graduate-level course based on the book "Predictability of Weather and Climate"{cite:p}`palmer2006predictability` will help you walk through some fundamental ideas in chaos (a.k.a. Butterfly effect), predictability and information theory. The goal is to learn how to (1) qualitatively and quantatively describe "predictability", (2) physically interpret dynamical systems and (3) find the optimal patterns favoring model error growth in various systems. Prerequisites for this course including: applied math (ODE, PDE and linear algebra), numerical analysis (Python or Julia), statistics and atmospheric dynamics.   


:::{note}
This is the v1.0 handout of the course "Chaos and Predictability" taught in the [Department of Atmospheric Science](http://www.as.ntu.edu.tw/index.php/eng), [National Taiwan University](https://www.ntu.edu.tw/chinese2007/english/index.html). There might be tons of typos. If you find any error, feel free to DM me (email Dr. Kai-Chih Tseng: kaichiht@princeton.edu)
:::

## Course Outlines
__Part I: Lecturing__
* {ref}`week1`
	* Climatological distribution vs forecast distribution 
    * State-dependent predictability and the mathematical assumptions 
    * Where the uncertainty comes from?
    * Introduction to ensemble forecast in weather and climate
    * HW1: Lorenz 63 model (testing state-dependent predictability) (10%)    
* {ref}`week2`
	* The origins of predictability (Scale Separation) and balanced dynamics
	* Atmospheric Blocking
	* Teleconnections 
	* HW2: Baroclinic instability and Predictability (15%)
* Week 3: __The minimalist models for studying chaos and predictability I: Lorenz 63__  
  	* The physical and mathematical background
    * Perfect model experitment
* Week 4: __The minimalist models for studying chaos and predictability II: Lorenz 96__
	* The physical and mathematical background
	* Single Scale
	* Multi-scale interations
	* Layapunov Exponent
	* The early and late stages of error growth in a multi-scale flow. 
	* HW3: Lorenz 96 model (calculating the Layapynov exponent) (15%)
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
	* HW4: Derivation of Liouville equation + 1st-order ODE case simulation (15%)
* Week 7: __The stability theory__
	* von Neumann method
	* Energy Method
	* Stability in Geophysical Fluid Dynamics (Shear instability and Eady problem)
* Week 8: __The generalized stability theory__
	* Generalized stability theorem
	* Modal vs non-modal systems
    * Reduced-dimension problems
* Week 9: __Information theory__
	* The connection between Layapunov exponent and perfect model experiments
	* Average Predictability Time
	* Predictable Components
	* HW5: Identifying the predictable components in full-physics operational forecast systems (20%)
* Week 10: __Future projection problems__
	* Initial value problems vs boundary value problems
	* Generalized stability and Predictable Components in a non-linear/boundary value problem
	* Hasselmann (1993) Optimal fingerprints for the detection of time-dependent climate change
	* Thompson et al. (2015) Quantifying the Role of Internal Climate Variability in Future Climate Trends




__Part II: Literature Reviews__ (35%)
* Week 11: __The predictability of a flow which possesses many scales of motions (Lorenz 1969)__
	* Key Lecture: Week 8	
* Week 12: __Stochastic Forcing of ENSO by the Intraseasonal Oscillation (Moore and Kleeman 1999)__
	* Key Lecture: Week 9
* Week 13: __Ensemble-based sensitivity analysis (Hakim and Torn 2008)__
	* Key Lecture: Week 9
* Week 14: __The Critical Role of Non-Normality in Partitioning Tropical and Extratropical Contributions to PNA Growth (Henderson et al. 2020)__
	* Key Lecture: Week 9
* Week 15: __Physically Interpretable Neural Networks for the Geosciences: Applications to Earth System Variability (Toms et al. 2020)__
	* Key Lecture: Week 10

## Grading
HW X 5 (10%, 15%, 15%, 15% 20%)
Literature Review (25%)


## Book
```{bibliography} references.bib
:filter: docname in docnames
```
