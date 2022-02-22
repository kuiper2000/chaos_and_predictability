# Week 1: Predictability of Weather and Climate

In this section, we will walk through a few fundamental concepts of predictability, including how does it connect to statistics, differential equations and linear algebra. You will find out the definition of predictability is relatively straightforward. 

__Predictability__
\
When we say something (e.g., x) is predictable, that means we can somewhat track its time evolution (Let's say it is $x(t)$). However, with the increase of forecast lead time (defined as the time difference between initial state and the final state), our confidence about what $x(t)$ should be will also decrease. Thus, it is physically more reasonable to use a probability density function (PDF $p(x(t))$) to describe $x(t)$ since we are talking about how "confident" we are. 
\
\
Now, let's imagine two very extreme cases. The first case is $t\rightarrow 0$. In this case, we are very certain what $x(t)$ should be. Thus, it's not hard to find $p(x(t))$ is very narrow and almost like a delta function (Fig. 1, solid line). The second case is that $t\rightarrow \infty$. In this case, we don't have any confidence. Thus the best way to estimate $x(t)$ is by randomly sampling from its historical values (aka "population") (Fig. 1, dashed line). Thus, from $t=0$ to $t\rightarrow 0$, we can find $p(x(t))$ is chaging its shape continuously from one PDF to the other.
\
\
   

