import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("seaborn-bright")                         # Setting the styles for graphing.
plt.rcParams["font.size"] = 12

m, g, d, delta, r, R, L0, L1, alpha, c, k, b, phi = sym.symbols('m, g, d, delta, r, R, L0, L1, alpha, c, k, b, phi')

z1e = sym.symbols('z1e')

m_value = 0.425
g_value = 9.81
d_value = 0.42
delta_value = 0.65
#r_value = 0.125
R_value = 53
#L0_value = 0.120
#L1_value = 0.025
#alpha_value = 1.2
c_value = 6815
k_value = 1880
#b_value = 10.4
phi_value = (42/180 * np.pi)


z3e = (k_value*(z1e-d_value) - m_value*g_value*np.sin(phi_value))**0.5 * (delta_value-z1e)**2
z3e /= c_value**0.5

xmin = d_value + (m_value*g_value*np.sin(phi_value)/k_value)

span = np.linspace(xmin, delta_value, 500)

results = []
for x in span:
    results.append(z3e.subs(z1e, x) * R_value)

index = results.index(max(results))

plt.plot(span, results)
plt.title("Question B1 - $V^e$ dynamics")
plt.xlabel('Equilibrium position of ball, $x^e$  (m)')
plt.ylabel('Equilibrium Voltage, $V^e$  ($V$)')
plt.axhline(y=results[index], color='black', linestyle='dashed')
plt.axvline(x=span[index], color='black', linestyle='dashed')
plt.annotate(("$x^e$ = %.3fm" % span[index] + "\n$V^e$ = 0.2V"), xy=(0.48, 0.012),
             bbox=dict(boxstyle="round", fc="0.9"))
plt.grid()
plt.show()
