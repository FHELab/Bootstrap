import sympy
from sympy import symbols, expand, poly
import numpy as np

import sys 

sys.setrecursionlimit(65537) 


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


x = symbols('x')
coeff = None
pre_poly = 1
for i in range(-63, 64): # poly = (x-63)* ... * (x+63)
  pre_poly = pre_poly*(x+i)



for i in range(1, 64): # loop through 32768+1, 32768+63
  tmp_poly = pre_poly*(x-32768) 
  xj = i+32768
  tmp_factor = 1

  # calculate coefficient factor
  for f in range(-63, 64):
    tmp_factor = tmp_factor * (xj-f)
    tmp_factor = tmp_factor % 65537
  for f in range(32768-63, 32768+64):
    if f != xj:
      tmp_factor = tmp_factor * (xj-f)
      tmp_factor = tmp_factor % 65537
  tmp_factor = modInverse(tmp_factor, 65537)


  # expand the poly
  for j in range(1, 64):
    if (i != j):
      tmp_poly = tmp_poly*(x-(32768-j))*(x-(32768+j))
    else:
      tmp_poly = tmp_poly*(x-(32768-j))

  p = poly(expand(tmp_poly))
  tmp_coeff = p.all_coeffs()
  tmp_coeff_lagrange = [(i * tmp_factor * 32768) % 65537 for i in tmp_coeff]
  if coeff is None:
    coeff = tmp_coeff_lagrange
  else:
    coeff = np.sum([coeff, tmp_coeff_lagrange], axis=0)
  coeff = np.mod(coeff, 65537)



for i in range(-63, 0): # loop through 32768-63, 32768-1
  tmp_poly = pre_poly*(x-32768) 
  xj = i+32768
  tmp_factor = 1

  # calculate coefficient factor
  for f in range(-63, 64):
    tmp_factor = tmp_factor * (xj-f)
    tmp_factor = tmp_factor % 65537
  for f in range(32768-63, 32768+64):
    if f != xj:
      tmp_factor = tmp_factor * (xj-f)
      tmp_factor = tmp_factor % 65537
  tmp_factor = modInverse(tmp_factor, 65537)


  # expand the poly
  for j in range(1, 64):
    if (i != -j):
      tmp_poly = tmp_poly*(x-(32768-j))*(x-(32768+j))
    else:
      tmp_poly = tmp_poly*(x-(32768+j))

  
  p = poly(expand(tmp_poly))
  print((p(xj) * tmp_factor * 32768) % 65537)
  tmp_coeff = p.all_coeffs()
  tmp_coeff_lagrange = [(i * tmp_factor * 32768) % 65537 for i in tmp_coeff]
  if coeff is None:
    coeff = tmp_coeff_lagrange
  else:
    coeff = np.sum([coeff, tmp_coeff_lagrange], axis=0)
  coeff = np.mod(coeff, 65537)



# add poly for x = 32768
tmp_factor = 1
tmp_poly = pre_poly
xj = 32768
for f in range(-63, 64):
  tmp_factor = tmp_factor * (xj-f)
  tmp_factor = tmp_factor % 65537
for f in range(32768-63, 32768+64):
  if f != xj:
    tmp_factor = tmp_factor * (xj-f)
    tmp_factor = tmp_factor % 65537
tmp_factor = modInverse(tmp_factor, 65537)


# expand the poly
for j in range(1, 64):
  tmp_poly = tmp_poly*(x-(32768-j))*(x-(32768+j))


p = poly(expand(tmp_poly))
tmp_coeff = p.all_coeffs()
print((p(xj) * tmp_factor * 32768) % 65537)
tmp_coeff_lagrange = [(i * tmp_factor * 32768) % 65537 for i in tmp_coeff]
if coeff is None:
  coeff = tmp_coeff_lagrange
else:
  coeff = np.sum([coeff, tmp_coeff_lagrange], axis=0)
coeff = np.mod(coeff, 65537)


print(coeff)


for i in range(32768-63, 32768+64):
  print(eval_poly(coeff[::-1], i) % 65537)

for i in range(-63, 64):
  print(eval_poly(coeff[::-1], i) % 65537)

