import scipy.optimize as opt
import matplotlib.pyplot as plt
import numpy as np

def f(x):
  return pow(2, x) - 3*x - 2

def fi(x):
  return np.log2(3*x + 2)
 
 
def Dichotomy(f, a, b, eps=1e-4):
  if f(a) * f(b) >= 0:
    return
  iter = 0
  print("Dichotomy method")
  while (abs(b - a) >= eps):
    x = (a + b) / 2.0
    if (f(a) * f(x) <= 0):
      b = x
    else:
      a = x
    iter += 1
  print(iter, "iterations")
  print("x = %0.5f" % x)

def Chord(f, a, b, eps=1e-4):
  if f(a) * f(b) >= 0:
    return
  iter = 1
  print("\nChord method")
  x0 = a
  xn = b - (f(b) / (f(b) - f(x0))) * (b - x0)
  while abs(xn - x0) > eps:
    x0 = xn
    xn = b - (f(b) / (f(b) - f(x0))) * (b - x0)
    iter += 1
  print(iter, "iterations")
  print("x = %0.5f" % xn)

def Iteration(f, a, b, eps=1e-4):
  iter = 0
  print("\nSimple iterations method")
  delta = 1
  x0 = b
  while delta > eps:
    if f(a)*f(b) >= 0:
      return 
    x = fi(x0)
    delta = abs(x-x0)
    x0 = x
    iter += 1
  print(iter, "iterations")
  print("x = %0.5f" % x)

Dichotomy(f, 2, 4)
Chord(f, 2, 4)
Iteration(f, 2, 4)

x = opt.fsolve(f, 4)[0]
print("\nSciPy\nx = %0.5f" % x)

x = np.arange(2, 4, 0.01)

plt.grid(True)
plt.plot(x, f(x), lw=2, color="blue")
plt.show()