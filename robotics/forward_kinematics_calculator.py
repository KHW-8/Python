import sympy as sp

alpha, a, d, theta = sp.symbols("α, a, d, Θ")

alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = sp.symbols("α1 α2 α3 α4 α5 α6")
a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
d1, d2, d3, d4, d5, d6 = sp.symbols("d1, d2, d3, d4, d5, d6")
theta1, theta2, theta3, theta4, theta5, theta6 = sp.symbols("Θ1 Θ2 Θ3 Θ4 Θ5 Θ6")

T01, T12, T23, T34, T45, T56 = sp.symbols("T01 T12 T23 T34 T45 T56");

R_X_alpha = sp.Matrix([
        [1,     0,              0,              0],
        [0,     sp.cos(alpha),  -sp.sin(alpha), 0],
        [0,     sp.sin(alpha),  sp.cos(alpha), 0],
        [0,     0,              0,              1]
    ])

D_X_a= sp.Matrix([
    [1, 0, 0, a],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])


R_Z_theta = sp.Matrix([
    [sp.cos(theta),    -sp.sin(theta),    0,  0],
    [sp.sin(theta),    sp.cos(theta),     0,  0],
    [0,                 0,                1,  0],
    [0,                 0,                0,  1],
])

D_Z_d = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, d],
    [0, 0, 0, 1]
])

class DH_Param:
    def __init__(self, alpha, a, d, theta):
        self.alpha = alpha
        self.a = a
        self.d = d
        self.theta = theta

def calculate(dh_params):
    list = []

    for dh_param in dh_params:
        R_X_alpha_tmp = R_X_alpha
        D_X_a_tmp = D_X_a
        R_Z_theta_tmp = R_Z_theta
        D_Z_d_tmp = D_Z_d

        for name, value in vars(dh_param).items():
            if isinstance(value, sp.Symbol):
                match name:
                    case "alpha":
                        R_X_alpha_tmp = R_X_alpha.subs(alpha, dh_param.alpha)
                    case "a":
                        D_X_a_tmp = D_X_a.subs(a, dh_param.a)
                    case "d":
                        D_Z_d_tmp = D_Z_d.subs(d, dh_param.d)
                    case "theta":
                        R_Z_theta_tmp = R_Z_theta.subs(theta, dh_param.theta)
            else:
                match name:
                    case "alpha":
                        R_X_alpha_tmp = R_X_alpha.evalf(chop=True, subs={alpha: dh_param.alpha})
                    case "a":
                        D_X_a_tmp = D_X_a.evalf(chop=True, subs={a: dh_param.a})
                    case "d":
                        D_Z_d_tmp = D_Z_d.evalf(chop=True, subs={d: dh_param.d})
                    case "theta":
                        R_Z_theta_tmp = R_Z_theta.evalf(chop=True, subs={theta: dh_param.theta})


        list.append(sp.nsimplify(R_X_alpha_tmp * D_X_a_tmp * R_Z_theta_tmp * D_Z_d_tmp))

    return list

def product(list, begin=0, end=0):
    M = sp.eye(4)

    if end != 0 and begin < end:
        for i in range(begin, end):
            M *= list[i]
    else:
        for T in list:
            M *= T

    M = sp.nsimplify(M)
    M = sp.expand(M)
    M = sp.trigsimp(M)

    return M

if __name__ == "__main__":

    list = calculate([
        DH_Param(0, 0, 0, theta1),
        DH_Param(0, a1, 0, theta2),
        DH_Param(0, a2, 0, theta3)
    ])

    for T in list:
        sp.pprint(T)

    T = product(list)
    sp.pprint(T)