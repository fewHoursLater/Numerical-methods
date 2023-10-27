It is necessary to numerically solve the integral equation with a given accuracy $epsilon$<br />
$u(x)-\int_{a}^{b}K(x,s)u(s)ds=f(x),  x\in [a,b]$<br />
Here $u(x)$ is the desired function with the domain of definition $[a,b]$. The given functions $K(x, s), u(x)$ are called, respectively, the kernel and the right side of the integral equation. <br />

I use elementary trapezoid formulas, as well as Gauss quadrature formulas for 3 nodes. I divide the initial segment into many small ones, on each of which I carry out numerical differentiation, the results of which are then summed up and the Runge rule is applied
