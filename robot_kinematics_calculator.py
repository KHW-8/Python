import sympy as sp

alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = sp.symbols("α1 α2 α3 α4 α5 α6")
a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
d1, d2, d3, d4, d5, d6 = sp.symbols("d1, d2, d3, d4, d5, d6")
theta1, theta2, theta3, theta4, theat5, theta6 = sp.symbols("Θ1 Θ2 Θ3 Θ4 Θ5 Θ6")

T01 = sp.Matrix([
    [sp.cos(), 0, 0, 0],
    [sp.sin(), 1, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]
])

T12 = sp.Matrix([
    []
])

sp.pprint(theta1)