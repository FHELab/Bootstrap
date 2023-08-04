import sympy
from sympy import symbols, expand, poly
import numpy as np

import sys 

prime = 65537
sys.setrecursionlimit(prime) 
x = symbols('x')

xs = []
ys = []

xo = [57004, 46969, 21931, 39030, 59092, 9965, 30013, 58301]
yo = [51237, 63408, 61771, 29236, 28266, 35197, 57032, 60306]
# xo = [57004]
# yo = [51237]

for i in range(len(xo)):
   for x_val in range(xo[i]-63, xo[i]+64):
      xs.append(x_val)
      ys.append(yo[i])

# for i in range(-63, 64):
#    xs.append(i)
#    ys.append(0)
# for i in range(32705, 32832):
#    xs.append(i)
#    ys.append(32768)


def power(x, y, m):
    if (y == 0):
        return 1
    p = power(x, y / 2, m) % m
    p = (p * p) % m

    if (y % 2 == 0):
      return p
    else:
      return (x * p) % m


def modInverse(a, m):
    return power(a, m - 2, m)


def eval_poly(lst, x): 
  n, tmp = 0, 0
  for a in lst:
    tmp = tmp + (a * (x**n))
    n += 1

  return tmp


def lagrange_interpolate(xs, ys):
    coeff = None
    for i in range(len(xs)):
        print(i)
        if ys[i] != 0:
            tmp_poly = 1
            tmp_factor = 1
            for x_value in xs:
                if x_value != xs[i]:
                    tmp_poly = tmp_poly * (x - x_value)
                    tmp_factor = (tmp_factor * (xs[i] - x_value)) % prime
            tmp_factor = modInverse(tmp_factor, prime)

            p = poly(expand(tmp_poly))
            tmp_coeff = p.all_coeffs()
            tmp_coeff_lagrange = [(c * tmp_factor * ys[i]) % prime for c in tmp_coeff]
            if coeff is None:
                coeff = tmp_coeff_lagrange
            else:
                coeff = np.sum([coeff, tmp_coeff_lagrange], axis=0)
            coeff = np.mod(coeff, 65537)

    return coeff




coeff = lagrange_interpolate(xs, ys)


for i in range(len(xs)):
  print(xs[i], eval_poly(coeff[::-1], xs[i]) % 65537, ys[i])


print(coeff)


# print(eval_poly([1,1,2], 0))
# print(eval_poly([1,1,2], 2))
