""" Here we compare our solver with published Schwarzschild QNM values (units
of 2M=1)

The known values taken from Berdi, Cardoso and Starinets "Quasinormal Modes
of Black Holes and Black Branes", Class. Quantum Grav. 16 163001 (2009),
arXiv: 0905.2975, as well as references mentioned in this paper :

scalar field l=0, s=0, n=0: M*omega ~ 0.1105 - 0.1049i , omega ~ 0.2210 -
0.2098i
scalar field l=1, s=0, n=0: M*omega ~0.2929 - 0.0977i , omega ~ 0.5858 -
0.1954i
gravitational field l=2, s=2, n=0: M*omega ~0.3737 - 0.0890i , omega ~0.7474 -
0.1779i
gravitational field l=2, s=2, n=1: M*omega ~0.3467 - 0.2739i , omega ~0.6934 -
0.5478i

Errors will be taken into account up to 1e-3, because we are mainly interested
in only large errors. The solver already has high accuracy and precision."""

import mpmath as mp

""" (l, s, n, guess, expected_omega) """
known_vals = [
 (0,0,0, mp.mpc(0.22, -0.21), mp.mpc(0.2210, -0.2098)),
 (1,0,0, mp.mpc(0.59, -0.20), mp.mpc(0.5858, -0.1954)),
 (2,2,0, mp.mpc(0.75, -0.18), mp.mpc(0.7474, -0.1779)),
 (2,2,1, mp.mpc(0.69, -0.55), mp.mpc(0.6934, -0.5478)),
 ]

def test():
  mp.mp.dps = 40
  for l, s, n, guess, expected_omega in known_vals:
    omega = qnms(l, s, n, guess, N=80)
    difference = abs(omega-expected_omega)
    print(f"l={l}, s={s}, n={n}: got {omega}, expected {expected_omega}, |difference|={difference}")
    assert difference < mp.mpf("1e-3"), f"QNM mismatch for l={l}, s={s}, n={n}"
if __name__ == "__main__":
   test()
   print("All known value tests are correct.")
