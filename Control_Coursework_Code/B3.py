import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

m = 0.425
g = 9.81
d = 0.42
delta = 0.65
r = 0.125
R = 53
L0 = 0.120
L1 = 0.025
alpha = 1.2
c = 6815
k = 1880
b = 10.4
phi = (42/180 * np.pi)

z1e = 0.467
ve = 0.2
z3e = 0.017

a1 = (2 * z3e) / ((delta - z1e)**2)
a2 = (2 * (z3e**2)) / ((delta - z1e)**3)
a3 = 1 / (L0+L1*np.exp(-alpha*(delta-z1e)))

tf_num = [c*a1*a3]
tf_den = [(7*m/5), (b+((7*m/5)*a3*R*m)), (-k-c*a2+a3*R*b), (-a3*R*k - c*a2*a3*R)]

Gx = ctrl.TransferFunction(tf_num, tf_den)

Gx_imp, y_imp = ctrl.impulse_response(Gx)
Gx_step, y_step = ctrl.step_response(Gx)

plt.plot(Gx_imp, y_imp)
plt.plot(Gx_step, y_step)
plt.show()