import numpy as np





c_p = 1004.0
theta0 = 300.0
g0     = 10.0
sigma  = 5.67e-8
H      = 1e4
rho0   = 1.0

B0 = 10.0

t = np.linspace(0, 1e2*86400.0, 1000)

C0 = (1 + B0 / g0)**(-3)
tmp = 3 * sigma * theta0**3 / (c_p * rho0 * H)


B = g0 * ( 1/( C0 + tmp * t )**(1/3) - 1)

print("Loading matplotlib")
import matplotlib.pyplot as plt
print("Done")

fig, ax = plt.subplots(1, 1)

ax.plot(t / 86400, B)

plt.show()


