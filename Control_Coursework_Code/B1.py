import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("seaborn-bright")                         # Setting the styles for graphing.
plt.rcParams["font.size"] = 12

z1e = sym.symbols('z1e')

m = 0.425
g = 9.81
d = 0.42
delta = 0.65
R = 53
c = 6815
k = 1880
phi = (42 / 180 * np.pi)

z3e = (k * (z1e - d) - m * g * np.sin(phi)) ** 0.5 * (delta - z1e) ** 2
z3e /= c ** 0.5

xmin = d + (m * g * np.sin(phi) / k)
span = np.linspace(xmin, delta, 500)

results = []
for x in span:
    results.append(z3e.subs(z1e, x) * R)

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
