import numpy as np
import matplotlib.pyplot as plt
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


def non_linear_system(state,t):    # Non-Linear System
    x = state[0]
    xd = state[1]
    I = state[2]

    v = R*(delta - x) * ((k * (x - d) - m * g * np.sin(phi)) ** 0.5)
    v /= (c ** 0.5)

    xdd = (5/7*m)*(m*g*np.sin(phi)+c*((I/(delta-x))**2) - k*(x-d)-b*xd)
    Id = (v - I*R) / (L0 + L1*np.exp(-alpha*(delta-x)))

    return [xd, xdd, Id]


x = 0.75*(d + ((m*g*np.sin(phi))/k)) + 0.25*delta

state0 = [x, 0, 0]
t = np.linspace(0.0, 1.0, 100)

state = odeint(non_linear_system, state0, t)

plt.plot(t, state)
plt.title("Question B2 - Non-Linear dynamics")
plt.xlabel('Time, t (s)')
plt.ylabel('STATES')
plt.legend(('$x$ $(m)$ ', '$\dot{x}$ $(ms^{-1})$', '$I$ $(A)$'))
plt.grid()
plt.show()
