# Week 1: Predictability of Weather and Climate

In this section, we will walk through a few fundamental concepts of predictability, including how does it connect to statistics, differential equations and linear algebra. You will find out the definition of predictability is relatively straightforward. 

__Predictability__
\
When we say something (let's say $x(t)$ in our case) is predictable, that means we can somewhat track its time evolution. However, in most cases, our confidence about what $x(t)$ should be will decrease with the increase of forecast lead time (defined as the time difference between initial state and the final state). Thus, it is physically more reasonable to use a probability density function (PDF $p(x(t))$) to describe $x(t)$ since we are talking about how "confident" we are. 
\
\
Now, let's imagine two very extreme cases. The first case is $t\rightarrow 0$. In this case, we are very certain what $x(t)$ should be. Thus, it's not hard to find $p(x(t))$ is very narrow and almost like a delta function (Fig. 1, solid line). The second case is that $t\rightarrow \infty$. In this case, we don't have any confidence. Thus the best way to estimate $x(t)$ is by randomly sampling from its historical values (aka "population") (Fig. 1, dashed line). Thus, from $t=0$ to $t\rightarrow \infty$, we can find $p(x(t))$ is chaging its shape continuously from one PDF to the other.
\
\
![FIG1](FIG1.png)
\
\
Here we can formulte our mathematical definition of "time of predictability limit". It is defined as the "moment" that the following null hypothesis is rejected.


$$ 
H_0: p(x(t)) = p(x(\infty))
$$ (my_other_label)

Eq. 1 indicates that the moment we can no longer tell the difference between $p(x(t))$ 