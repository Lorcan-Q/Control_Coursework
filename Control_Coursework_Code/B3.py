import control as ctrl
import numpy as np
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

Gx = ctrl.series(a1_val*a4_val, tf1, tf2)

Gx_imp, y_imp = ctrl.impulse_response(Gx)
Gx_step, y_step = ctrl.step_response(Gx)

plt.plot(Gx_imp, y_imp)
plt.grid()
plt.title("Question B3 - Impulse Response")
plt.xlabel('Time, t (s)')
plt.ylabel('Equilibrium Displacement, $\overline{X}$ ($m$)')
plt.show()

plt.plot(Gx_step, y_step)
plt.title("Question B3 - Step Response")
plt.xlabel('Time, t (s)')
plt.ylabel('Equilibrium Displacement, $\overline{X}$ ($m$)')
plt.grid()
plt.show()