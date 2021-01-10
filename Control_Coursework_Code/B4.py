import control as ctrl
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

plt.style.use("seaborn-bright")                         # Setting the styles for graphing.
plt.rcParams["font.size"] = 12

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

xe = 0.75*(d + ((m*g*np.sin(phi))/k)) + 0.25*delta
ve = R*(delta - xe) * ((k * (xe - d) - m * g * np.sin(phi)) ** 0.5)
ve /= (c ** 0.5)
Ie = ve/R

a1_val = (10 * c * Ie) / (7 * m * ((delta - xe) ** 2))
a2_val = (5 / (7 * m)) * (-k + ((2 * c * (Ie ** 2)) / ((delta - xe) ** 3)))
a3_val = (5 * b) / (7 * m)
a4_val = 1 / (L0 + (L1 * np.exp((-alpha * (delta - xe)))))

tf1 = ctrl.TransferFunction(1, [1, (a4_val*R)])
tf2 = ctrl.TransferFunction(1, [1, a3_val, -a2_val])

Gx = a1_val * a4_val * ctrl.series(tf1, tf2)

ctrl.bode(Gx, dB=True)  # LFA appx = -10dB HFA appx.= -10 - 60log(w)
plt.suptitle("B4 - Bode Plot")
plt.show()
