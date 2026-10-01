import sympy as sp
import math

alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = sp.symbols("α1 α2 α3 α4 α5 α6")
a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
d1, d2, d3, d4, d5, d6 = sp.symbols("d1, d2, d3, d4, d5, d6")
theta1, theta2, theta3, theta4, theta5, theta6 = sp.symbols("Θ1 Θ2 Θ3 Θ4 Θ5 Θ6")

link1 = 0.1
link2 = 0.1
link3 = 0.1
link4 = 0.1


px = a2*sp.cos(theta1)*sp.cos(theta2) + a3*sp.cos(theta1)*sp.cos(theta2 + theta3) - d2*sp.sin(theta1) + d5*sp.sin(theta2 + theta3 + theta4)*sp.cos(theta1)
py = a2*sp.sin(theta1)*sp.cos(theta2) + a3*sp.sin(theta1)*sp.cos(theta2 + theta3) - d2*sp.cos(theta1) + d5*sp.sin(theta2 + theta3 + theta4)*sp.cos(theta1)

px = px.evalf(chop=True, subs={
    theta1: math.radians(10),
    a2: link2,
    a3: link3,
    d2: link1,
    d5: link4,
})

py = py.evalf(chop=True, subs={
    theta1: math.radians(10),
    a2: link2,
    a3: link3,
    d2: link1,
    d5: link4,
})

sp.pprint(px)
sp.pprint(py)

expr1 = sp.expand(px**2) 
expr2 = sp.expand(py**2)
sp.pprint(expr1)
sp.pprint(expr2)

# expr = sp.simplify(expr1 + expr2)
# sp.pprint(expr)