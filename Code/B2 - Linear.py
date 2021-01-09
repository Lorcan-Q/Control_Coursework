import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy.integrate import odeint

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


def linearSystem(state,t):    # Non-Linear System
    x = state[0]
    xd = state[1]
    I = state[2]

    v = R*(delta - x) * ((k * (x - d) - m * g * np.sin(phi)) ** 0.5)
    v /= (c ** 0.5)

    a1_val = (2 * I) / ((delta - x) ** 2)
    a2_val = (2 * (v/R ** 2)) / ((delta - x) ** 3)
    a3_val = 1 / (L0 + L1 * np.exp(-alpha * (delta - x)))

    xdd = (5/7*m)*(c*(a1_val*(I - v/R) + a2_val*x) - k*x - b*xd)
    Id = a3_val*(-R*I + v)

    return [xd, xdd, Id]


x = 0.75*(d + ((m*g*np.sin(phi))/k)) + 0.25*delta

state0 = [x, 0, 0]
t = np.linspace(0.0, 1.0, 100)

state = odeint(linearSystem, state0, t)

plt.plot(t, state)
plt.legend(('$x$', '$\dot{x}$', 'I'))
plt.grid()
plt.show()