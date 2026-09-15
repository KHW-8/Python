import sympy as sp

alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = sp.symbols("α1 α2 α3 α4 α5 α6")
a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
d1, d2, d3, d4, d5, d6 = sp.symbols("d1, d2, d3, d4, d5, d6")
theta1, theta2, theta3, theta4, theat5, theta6 = sp.symbols("Θ1 Θ2 Θ3 Θ4 Θ5 Θ6")

# T01
T01 = sp.Matrix([
    [sp.cos(theta1),     -sp.sin(theta1),   0,          0],
    [sp.sin(theta1),    sp.cos(theta1),     0,          0],
    [0,                 0,                  1,          0],
    [0,                 0,                  0,          1]
])

# T12
M1 = sp.Matrix([
    [1,     0,                  0,            ],
    [0,     sp.cos(alpha1),     sp.sin(alpha1)],
    [0,     -sp.sin(alpha1),    sp.cos(alpha1)]
])

M2 = sp.Matrix([
    [sp.cos(theta2),    -sp.sin(theta2),    0],
    [sp.sin(theta2),    sp.cos(theta2),     0],
    [0,                 0,                  1]
])

M1 = M1.evalf(chop=True, subs={alpha1: sp.rad(90)})

M3 = M1 * M2

M3 = M3.row_insert(3, sp.Matrix([[0, 0, 0]]))
T12 = M3.col_insert(3, sp.Matrix([0, d2, 0, 1]))


# T23


#
sp.pprint(sp.nsimplify(T01))