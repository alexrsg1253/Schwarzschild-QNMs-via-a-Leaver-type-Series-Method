"""
Schwarzschild QNMs via a Leaver-type Series Method

Setup:

-Units: r=2M=1 (Schwarzschild radius: event horizon)

-Regge-Wheeler equation for a spin-s perturbation:

[d^2/dr*^2 + w^2-V(r)] psi =0,

where V(r)=f(r)[l(l+1)/r^2 + (1-s^2)/r^3]   and  f(r)= 1-1/r .

Here, l is the angular momentum number and r* is the tortoise coordinate.
s=0 corresponds to a scalar field, s=1 to an electromagnetic field and s=2
to a gravitational field.

Now we define a variable z=1-1/r, such that z=0 corresponds to the event
horizon (r=1) and z=1 is spatial infinity.

By factoring out the two singular prefactors from the physical boundary
conditions, 

psi(r) = e^(i w r) * z^(-i w) * phi(z),

we get a Fuchsian ODE for phi with polynomial coefficients,

A(z) phi'' + B(z) phi' + C(z) phi = 0

A, B and C are polynomials of degree 4, 3 and 2, respectively, and their 
coefficients depend on w, l and s. 

Writting an ansatz for phi as a power series solution,

phi(z) = \sum_n a_n z^n

plugging it into the ODE and matching powers of z gives a linear recursion 
relating a_(n+1), a_n, a_(n-1) and a_(n-2). 

For n<0, a_n=0, so the physical solution corresponds to a normalization choice
of a_0=1. Also, the recessive solution of the recursion relation must be 
picked out, since it represents the boundary condition at infinity (purely 
outgoing waves). 

Following Leaver's approach, this condition is the same as requiring that a 
large truncated version of the recursion matrix be singular. To get the QNM 
frequencies, w, we find the roots on the determinant det(M_N(w)) for increasing
truncation size N and then checking for convergence.

This procedure will be validated using already published Schwarzschild QNM 
tables for various values of (l, s, n).
"""

import mpmath as mp

def poly_coeffs(w, l, s):
  """Polynomial coefficients A(z), B(z) and C(z) ordered in powers of z from
  lower to higher"""
  i = mp.mpc(0,1)
  A= [mp.mpc(0), mp.mpc(-1), mp.mpc(3), mp.mpc(-3), mp.mpc(1)]
  B= [2*i*w-1, -8*i*w+5, 8*i*w-7, -2*i*w+3]
  C= [l**2+l+1-s**2-6*w**2-3*i*w, -l**2-l-2+2*s**2+5*w**2+5*i*w, 
      1-s**2-w**2-2*i*w]
  return A, B, C

def recur_coeffs(w, l, s, n):
  """Coefficients of the recursion relation given by alpha, beta, gamma, 
  delta at index n: 
  alpha*a_(n+1) + beta*a_n + gamma*a_(n-1) + delta*a_(n-2)=0"""
  A, B, C = poly_coeffs(w, l, s)
  alpha = A[1]*(n+1)*n + B[0]*(n+1) 
  beta = A[2]*n*(n-1) +B[1]*n+C[0]
  gamma = A[3]*(n-1)*(n-2)+B[2]*(n-1)+C[1]
  delta = A[4]*(n-2)*(n-3)+ B[3]*(n-2)+C[2]
  return alpha, beta, gamma, delta

def matrix(w, l, s, N):
  """ NxN row-normalized truncated recursion matrix. A row n has,
  a_(n+1)+(beta/alpha)a_n+(gamma/alpha)a_(n-1)+(delta/alpha)_a(n-2)=0.
  We are row-normalizing with alpha, since it's the coefficient of the 
  leading term. It keeps the entries of the matrix of order 1 for large n,
  which ensures the determinant won't diverge"""
  M = mp.matrix(N,N)
  for n in range(N):
    alpha, beta, gamma, delta = recur_coeffs(w, l, s, n)
    if alpha == 0:
       alpha = mp.mpf("1e-30")
    row = {n+1: mp.mpc(1), n: beta/alpha, n-1: gamma/alpha, n-2: delta/alpha}
    for index, coeff in row.items():
      if 0 <= index < N:
        M[n, index] += coeff
  return M

  def hill_det(w, l, s, N):
    """The determinant of the matrix M, det(M_N(w)), vanishes in the limit 
    where N goes to infinity at the QNM frequencies"""
    return mp.det(matrix(mp.mpc(w), l, s, N))
  def qnms(l, s, n_overtone, omega, N=80, dps=50):
    """Gives a QNM frequency at a certain value of omega. 
    l is an integer and it's the angular momentum (>= s). 
    s is an integer, spin of the perturbing field (0 scalar field, 1 EM field, 
    2 gravitational field). 
    n_overtone is an integer, the overtone number (n=0 is fundamental overtone).
    omega is for an estimated complex value of omega. 
    N is truncation size of Hill determinant. 
    dsp is from mpmath for decimal place precision. 
    This function returns the converged complex frequency omega."""
    mp.mp.dps = dps
    def f(w):
      return hill_det(w, l, s, N)
    return mp.findroot(f, mp.mpc(omega))

    def converged_qnm(l, s, n_overtone, omega, N_list=(40, 60, 80, 100, 130),
                     dps=50):
     """Same as qnms but it uses increases truncation number N. Each result
    is used as the next initial value for omega and it returns the sequence
    of estimates so that the convergence can be directly checked."""
     mp.mp.dps = dps
     root = mp.mpc(omega)
     empty = []
     for N in N_list:
       def f(w, N=N):
         return hill_det(w, l, s, N)
       root = mp.findroot(f, root)
       empty.append((N, root))
     return empty
