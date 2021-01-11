import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy.integrate import odeint

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
phi = ((42/180) * np.pi)

xe = (0.75 * (d + ((m*g*np.sin(phi))/k))) + (0.25*delta)


def linear_system(state, t):
    x = state[0]
    xd = state[1]
    I = state[2]

    ve = R * (delta - xe) * ((k * (xe - d) - m * g * np.sin(phi)) ** 0.5)
    ve /= (c ** 0.5)
    Ie = ve/R

    a1_val = (10 * c * Ie) / (7 * m * ((delta-xe)**2))
    a2_val = (5/(7*m)) * (-k + ((2*c*(Ie**2))/((delta-xe)**3)))
    a3_val = (5*b)/(7*m)
    a4_val = 1 / (L0 + (L1 * np.exp((-alpha * (delta-xe)))))

    xdd = (a1_val*(I - Ie)) + (a2_val*(x-xe)) - (a3_val*(xd-0))
    Id = a4_val * ((ve-ve) - R * (I-Ie))

    return [xd, xdd, Id]


state0 = [xe, 0, 0]
t = np.linspace(0.0, 1.0, 400)

state = odeint(linear_system, state0, t)

plt.plot(t, state)
plt.title("Question B2 - Linear dynamics")
plt.xlabel('Time, t (s)')
plt.ylabel('STATES')
plt.legend(('$x$ $(m)$ ', '$\dot{x}$ $(ms^{-1})$', '$I$ $(A)$'))
plt.grid()
plt.show()