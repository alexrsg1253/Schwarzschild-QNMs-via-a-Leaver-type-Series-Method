# Schwarzschild QNMs via a Leaver-type Series Method
Here I will use a Frobenius-series (Leaver-type) solver for Schwarzschild black hole (BH) quasinormal modes (QNMs), which will be Step 1 of a project that aims to compare this numerical method to the semiclassical CFT/Nekrasov-Shatashvili quantization condition (this will form part of Step 2). 
## Python Scripts
- `QNM_Leaver/schwarzschild.py` : The solver. In this script, the Frobenius recursion relation is derived for the Regge-Wheeler equation with spin s=0,1,2 and angular momentum number l. A truncated, row-normalized matrix M is constructed and the QNM frequencies are found from the roots of its determinant. Specifically, a Hill determinant method was used, similar to Leaver's continued-fraction method.
- `tests/test_known_values.py` : Regression tests againts published Schwarzschild QNM values for scalar and gravitational perturbations with the fundamental mode and the first overtone considered. 
- `example_computation.py` : An example that shows how the solver works and how to check for convergence as the truncation size N increases. 
## Starting Point
```bash
pip install -r Requirements.txt
PYTHONPATH=. python3 example_computation.py
PYTHONPATH=. python3 tests/test_known_values.py
```

The expected output for the fundamental gravitational mode (for l=2 and n=0) is

``
M*omega = 0.373672 - 0.088962i
``
which matches with published values (see Berti, Cardoso and Starinets, Class. Quantum Grav. 26.163001 (2009), arXiv: 0905.2975) up to 5/6 significant figures. 

## Method

The Regge-Wheeler equation is transformed by defining the variable `z=1-1/r` (so `z=0` is the horizon and `z=1` is spatial infinity). The two singular prefactors that come from the physical boundary conditions (purely ingoing at the horizon and purely outgoing at spatial infinity) are factored out, 

`` psi(r) = exp(i*omega*r) * z**(-i*omega) * phi(z) ``

Plugging this into the ODE gives a Fuchsian equation for `phi(z)` with polynomial coefficients. Using the power series ansatz `phi(Z)=sum_n a_n z^n`and matching powers of `z` yields a linear recursion relation which connects the coefficients `a_(n+1), a_n, a_(n-1), a_(n-2)`. The solution is recessive, and thus it is picked out by requiring a large truncated version of the recursion matrix to be singular (according to Leaver's method). QNM frequencies are the values of omega for which this occurs, and they're found by root-finding the determinant of the matrix and checking for convergence as the truncation size N grows. 

We don't have a three-term recursion, as expected, we get a four-term recursion. Despite this, the physics is the same and so are the boundary conditions of the system, and it has been checked against known QNM values successfully. 
