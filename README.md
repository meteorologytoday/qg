# Moist two-layer shallow water model

A moist two-layer shallow water model written in Julia. Reference paper:

Bembenek, E., Straub, D. N., & Merlis, T. M. (2020). Effects of moisture in a two-layer model of the midlatitude jet stream. Journal of the Atmospheric Sciences, 77(1), 131-147.

# Dynamics

The model construction references:

- Bouchut, F., & Zeitlin, V. (2010). A robust well-balanced scheme for multi-layer shallow water equations. Discrete and Continuous Dynamical Systems-Series B, 13(4), 739-758.



# Todo

1. Env and Grid
2. AtmState
3. OcnState
4. Construct dynamic core matrix and matrix operator
5. Euler backward method / RK4
6. Forward stepping and output
