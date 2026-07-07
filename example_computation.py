""" This is an example for computing the gravitational Schwarzschild QNM 
corresponding to the fundamental mode. We also check for convergence as 
truncation size N is increased. 
We run this with PYTHONPATH=. python3 example_computation.py """

from qnm_leaver import converged_qnm

if __name__ == "__main__":
  l, s, n =2, 2, 0
  guess = 0.75-0.18j
  print(f"Computing Schwarzschild QNM: l={l}, s={s}, n={n}")
  print("Convergence as truncation size N increases:\n")

  empty= converged_qnm(l, s, n, guess)
  for N, omega in empty:
    print(f"N={N:4d}: omega={omega}")
  
  print("n\In standard M*omega units, this is:")
  N_final, omega_final = empty[-1]
  print(f"M*omega={omega_final /2}")
