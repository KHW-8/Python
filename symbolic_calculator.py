import sympy as sp

c1 = sp.cos(sp.symbols("1"))
c2 = sp.cos(sp.symbols("2"))
c23 = sp.cos(sp.symbols("23"))

s1 = sp.sin(sp.symbols("1"))
s2 = sp.sin(sp.symbols("2"))
s23 = sp.sin(sp.symbols("23"))

px, py = sp.symbols("px py")
a2, a3 = sp.symbols("a2 a3")
d3, d4 = sp.symbols("d3 d4")

exp = (-s1*px + c1*py)**2 + (c1*px + s1*py)**2 + (-px)**2

exp = sp.expand(exp)
exp = sp.nsimplify(exp)
exp = sp.trigsimp(exp)


sp.pprint(exp)

exp = d3**2 + (a3*c23 - d4*s23 + a2*c2)**2 + (a3*s23 + d4*c23 + a2*s2)**2

exp = sp.expand(exp)
exp = sp.trigsimp(exp)

sp.pprint(exp)