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

Gx = a1_val*a4_val*ctrl.series(tf1, tf2)    # System Transfer function

Gm = ctrl.TransferFunction(1, [0.03, 1])    # Sensor Transfer function

Kp = 1
Ki = 0
Kd = 0

def pid(kp, ki, kd):
    """
    This function constructs a PID controller; returning its associated transfer function, based upon
    the provided values for its gain components.
    :param kp: The proportional gain value.
    :param ki: The integral gain value.
    :param kd: The differential gain value.
    :return: The PID's transfer function based on the gain values provided.
    """
    diff = ctrl.TransferFunction([1, 0], 1)     # Defines the differential's transfer function.
    intgr = ctrl.TransferFunction(1, [1, 0])    # Defines the integral's transfer function.
    pid_tf = kp + (kd * diff) + (ki * intgr)    # Determines the PID's overall transfer function.
    return pid_tf


G_controller = -pid(Kp, Ki, Kd)  # The PID controller's transfer function.

closed_loop_tf = ctrl.feedback(ctrl.series(Gx, Gm), G_controller)

t_final = 1                                     # The total time the simulation will execute.
num_points = 1000                               # The number of sampling points that will be taken.
t_span = np.linspace(0, t_final, num_points)    # The array of times at which the samples will be taken.

x_res, t, _ = ctrl.forced_response(Gx, t_span, 1)

plt.plot(x_res, t)
plt.grid()
plt.show()