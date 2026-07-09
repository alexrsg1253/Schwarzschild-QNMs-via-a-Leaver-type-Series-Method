""" CFT/Nekrasov-Shatashvili quantization condition for QNMs with the Leaver-type
solver we made in qnm_leaver/schwarzschild.py.

The goal is to use the CFT/ instanton counting quantization condition for QNMs
and check it against the numerical method involving the Leaver-type values
computed here.

Key references that I use in this project are the papers:

Bonelli, Iossa, Panea Lichtig, Tanzini, "Exact solution of Kerr black hole
perturbations via CFT2 and instanton counting: Greybody factor,quasinormal
modes and Love numbers", Phys. Rev. D 105, 044047 (2022), arXiv:2105.04483.
In this paper the QNM quantization condition is described and applied to Kerr,
along with greybody factors and Love numbers.

Bonelli, Iossa, Panea Lichtig, Tanzini, "Irregular Liouville Correlators and
Connection Formulae for Heun Functions", Commun. Math. Phys. 397, 635 (2023),
arXiv:2201.04491. In this paper, the connection formulae for the Heun functions
and its confluences are derived explicitly.

Aminov, Grassi, Hatsuda, "Black Hole Quasinormal Modes and Seiberg-Witten
Theory", Annales Henri Poincare 23, 1951 (2022), arXiv:2006.06111. Here, the
same type of quantization condition is derived as in the first paper, and I use
it as a cross-check for conventions.

Lisovyy & Naidiuk, "Perturbative connection formulas for Heun equations",
J. Phys. A 55, 434005 (2022), arXiv:2208.01604. Finally, this paper is useful
when it comes to the perturbative expansion of the connection coefficients and
accessory parameter, which is what actually will be coded as a quantization
condition. 

Dictionary: 

lambda3 = -16i M w
E = -l(l+1) +8M^2 w^2-1/4
m1 = s-2i M w, m2= -s -2i M w, m3= -2iMw

Quantization condition:

pi_B^3 (E, m , lambda3, hbar=1) = 2pi(n+1/2), n=0,1,2,...

Matone relation for Nf=3:

E = a^2 - lambda3*dF_inst/dlambda3

Full quantum period for Nf=3:

pi_B= dF/da = -2alog(lambda3/4 - pi)-2ilog(Gamma(1+2ia)/Gamma(1-2ia))
-i*sum_j log(Gamma(1/2+mj-ia)/Gamma(1/2+mj+ia)) + dF_inst/da

Nekrasov instanton partition function: Combinatorial sum over pairs of 
Young tableaux. 
F_inst = -lim_(eps2 -> 0) eps2*logZ(ia, m, hbar=1, eps2)

Here F_inst is found through a numeric Richardson extrapolation.


Status:

-The combinatorial elements like arm/leg length and guage/matter factors
have been validated by checking that the limit where eps2 goes to 0 is finite 
at 1 and 2 instanton order. 
-The numeric Richardson extrapolation has been validated by comparing it with
the symbolic result at 1 and 2 instanton order and it agreed up to 10+ 
significant figures. 
-The quantization process runs successfully. 
-It has been tested against the l=2, s=2, n=0 fundamental mode, and the 
residual goes from around 12.6 to around 2.7, so it's not converged yet.
At this point (|lambda3/4|~1.5), the source paper needed up to 12 instanton
orders to get good agreement in other regimes. So, I'm working on that. 


"""

import mpmath as mp

mp.mp.dps = 30
N_max_def = 7
Eps2_list_def = [mp.mpf('0.1')/2**k for k in range(5)]

#Young Tableaux combinatorics

def partitions(n):
  """All integer partitions of n as non-increasing tuples"""
  if n == 0:
    return[()]
  result = []

  def helper(remaining, max_part, current):
    if remaining == 0:
      result.append(tuple(current))
      return
    for k in range(min(remaining, max_part), 0, -1):
      helper(remaining -k, k, current +[k])

  helper(n, n, [])
  return result 

def conjugate(Y):
  if not Y:
    return ()
  m=Y[0]
  return tuple(sum(1 for part in Y if part >= j) for j in range(1, m+1))

def boxes(Y):
  return [(i,j) for i, row in enumerate(Y, start=1) for j in range(1, row+1)]


def arm(Y, s):
  """Boxes to the right of s in Y"""
  i, j = s
  Yi = Y[i-1] if i <= len(Y) else 0
  return Yi - j

def leg(Y, s):
  """Boxes below s in Y"""
  i, j = s
  Yc = conjugate(Y)
  Ycj = Yc[j-1] if j <= len(Yc) else 0
  return Ycj -1

_pair_cache = {}

def pairs(n):
  """Pairs (Y1, Y2) of partitions with |Y1|+|Y2| = n. It's cached since it's 
  independent of a, m, epsilon1 and epsilon2"""
  if n not in _pair_cache:
    _pair_cache[n] = [
        (Y1, Y2) for n1 in range(n+1)
        for Y1 in partitions(n1) for Y2 in partitions(n-n1)
    ]
  return _pair_cache[n]

#Nekrasov gauge/matter factors

def z_gauge(Y1, Y2, a, e1, e2):
  alpha = {1: a, 2: -a}
  Ys = {1: Y1, 2: Y2}
  result = mp.mpc(1)
  for I in (1, 2):
    for J in (1, 2):
      YI, YJ = Ys[I], Ys[J]
      for s in boxes(YI):
        result /= (alpha[I]-alpha[J]-e1*leg(YJ, s)+e2*(arm(YI,s)+1))
      for s in boxes(YJ):
        result /= (alpha[I]-alpha[J]+e1*(leg(YI, s)+1)-e2*(arm(YJ,s)))

  return result

def z_matter(Y1, Y2, a, m_list, e1, e2):
  alpha = {1: a, 2: -a}
  Ys = {1: Y1, 2: Y2}
  result = mp.mpc(1)
  for mk in m_list:
    for I in (1,2):
      for (i,j) in boxes(Ys[I]):
        result *= (alpha[I]+mk+(i-mp.mpf(0.5))*e1 + (j - mp.mpf(0.5))*e2)
  return result

def z_tilde(N_max, a, m_list, e1, e2):
  """z_tilde_n is a coefficient of q^n in Z bedore the Lambda^(4-Nf/4)^n 
  prefactor for n=1...N_max"""
  out = {}
  for n in range(1, N_max +1):
    total = mp.mpc(0)
    for Y1, Y2 in pairs(n):
      total += z_gauge(Y1, Y2, a, e1, e2)*z_matter(Y1, Y2, a, m_list, e1, e2)
    out[n] = total
  return out 

def cumulants(N_max, Zt):
  """L_n (coefficients of logZ) from z_tilde_1... z_tilde_N through the log of
  power series recursion: L_m = Zt_m - (1/m) sum_(k<m)kL_k Zt_(m-k)"""
  L = {}
  for m_ in range(1, N_max +1):
    s = Zt[m_]
    for k in range(1, m_):
      s -= (mp.mpf(k) / m_)*L[k]*Zt[m_ - k]
    L[m_] = s
  return L

def f_richardson(N_max, a_phys, hbar, m_list, eps2_list=None):
  """f_n(a, hbar, m) = -hbar * lim_(eps2 to 0)eps2 * L_n(i*a, m, hbar, eps2) 
  for n=1... N_max through Richardson extrapolation. f_n are the Taylor 
  coefficients of F_inst in powers of q=Lambda_Nf/4"""
  if eps2_list is None:
    eps2_list = Eps2_list_def
  a = 1j * a_phys
  seqs = {n: [] for n in range(1, N_max+1)}
  for e2 in eps2_list:
    Zt = z_tilde(N_max, a, m_list, hbar, e2)
    L = cumulants(N_max, Zt)
    for n in range(1, N_max + 1):
      seqs[n].append(e2*L[n])
  f = {}
  N = len(eps2_list)
  for n in range (1, N_max +1):
    T = [seqs[n][:]]
    for level in range(1, N):
      prev = T[level-1]
      T.append([(2**level*prev[i+1]-prev[i])/(2**level -1) for i in range(N
      - level)])                                                                    
    f[n] = -hbar * T[N-1][0]
  return f


# F_inst and its derivatives

def F_inst(f, Lambda_Nf, N_max, use_pade=True):
  """F_inst(Lambda)= sum_n f_n q^n , q=Lambda/4"""
  q = Lambda_Nf / 4
  G_coeffs = [f[n] for n in range(1, N_max + 1)]
  if use_pade:
    L = (N_max - 1) // 2
    M = (N_max -1) - L
    try:
      p, qden = mp.pade(G_coeffs, L, M)
      num = sum(p[i]* q**i for i in range(len(p)))
      den = sum(qden[i]* q**i for i in range(len(qden)))
      return q * num / den
    except Exception:
      pass
  return q*sum(G_coeffs[i]* q**i for i in range(len(G_coeffs)))


def dFinst_dLambda(f, Lambda_Nf, N_max):
  h = mp.mpf('1e-12') * (abs(Lambda_Nf) + 1)
  return (F_inst(f, Lambda_Nf + h, N_max) - F_inst(f, Lambda_Nf - h, N_max))/(2*h)


def dF_da(a_phys, Lambda_Nf, hbar, m_list, N_max=N_max_def, eps2_list=None):
  """full quantum period Pi_B = dF/da, for Nf=3"""
  log_term = -2*a_phys*mp.log(Lambda_Nf/4)-hbar*mp.pi
  log_term += -2j * hbar * mp.log(mp.gamma(1+2j*a_phys/hbar)/mp.gamma(1-
  2j*a_phys/hbar))
  ferm = mp.mpc(0)
  for mj in m_list:
    ferm += mp.log(mp.gamma(mp.mpf(0.5)+(mj-1j*a_phys)/hbar)
     /mp.gamma(mp.mpf(0.5)+(mj+1j*a_phys)/hbar))
  log_term += -1j * hbar * ferm

  h = mp.mpf('1e-7')
  f_p = f_richardson(N_max, a_phys+h, hbar, m_list, eps2_list)
  f_m = f_richardson(N_max, a_phys+h, hbar, m_list, eps2_list)
  dFinst_da = (F_inst(f_p, Lambda_Nf, N_max)-F_inst(f_m, Lambda_Nf, N_max))
  / (2*h)
  return log_term + dFinst_da

#Dictionary and quantization condition

def dictionary(M, l, s, w):
  """Schwarzschild parameters to gauge theory parameters"""
  Lambda3 = -16j * M * w
  E = -l * (l+1) + 8*M**2 * w**2 - mp.mpf(0.25)
  m = [s - 2j*M*w, -s-2j*M*w, -2j*M*w]
  return E, m, Lambda3

def solve_a(E, Lambda3, hbar, m_list, a_guess, N_max=N_max_def, eps2_list=None):
  """Invert the equation of E to get a"""
  def eq(a_phys):
    f = f_richardson(N_max, a_phys, hbar, m_list, eps2_list)
    dFdL = dFinst_dLambda(f, Lambda3, N_max)
    return a_phys**2 - Lambda3 * dFdL - E
  return mp.findroot(eq, a_guess)


def quantization_res(w, M, l, s, n_overtone, a_guess, N_max=N_max_def, 
eps2_list=None):
  """Residual of the quantization condition. It should be 0 at true QNM
  frequencies (once N_max is big enough for the given prefactor |Lambda3/4|)"""
  E, m, Lambda3 = dictionary(M, l, s, w)
  a0 = solve_a(E, Lambda3, 1, m, a_guess, N_max, eps2_list)
  PiB = dF_da(a0, Lambda3, 1, m, N_max, eps2_list)
  target = 2 * mp.pi * (n_overtone + mp.mpf(0.5))
  return PiB - target

if __name__ == "__main__":
  """Sanity check against l=2, s=2, n=0 mode"""
  M = mp.mpf(0.5)
  l, s, n_overtone = 2, 2, 0
  w_known = mp.mpc(0.7473, -0.1779)
  residual = quantization_res(w_known, M, l, s, n_overtone, 
  a_guess=mp.mpc(1.0, 0.0))
  print(f"Residual at known Leaver frequency: {residual} (|.|={abs(residual)})")
  print("See module docstring 'Status' for what this demonstrates.")
