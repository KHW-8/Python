import sympy as sp
import math

alpha, a, d, theta = sp.symbols("α, a, d, Θ")

alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = sp.symbols("α1 α2 α3 α4 α5 α6")
a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
d1, d2, d3, d4, d5, d6 = sp.symbols("d1, d2, d3, d4, d5, d6")
theta1, theta2, theta3, theta4, theta5, theta6 = sp.symbols("Θ1 Θ2 Θ3 Θ4 Θ5 Θ6")

T01, T12, T23, T34, T45, T56 = sp.symbols("T01 T12 T23 T34 T45 T56");

R_X_alpha = sp.Matrix([
        [1,     0,              0,              0],
        [0,     sp.cos(alpha),  -sp.sin(alpha), 0],
        [0,     sp.sin(alpha),  sp.cos(alpha),  0],
        [0,     0,              0,              1]
    ])

D_X_a = sp.Matrix([
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
        for i in range(begin, end + 1):
            M *= list[i]
    else:
        for T in list:
            M *= T

    M = sp.nsimplify(M)
    M = sp.trigsimp(M)

    return M


def pulse_to_deg(pulse):
    return pulse * 240 / 1000

def inverse_transformation(T_B2A):
    R_B2A = T_B2A[0:3, 0:3]
    P_B2A = T_B2A[0:3, 3:4]

    R_A2B = sp.transpose(R_B2A)
    P_A2B = -R_A2B * P_B2A

    T_A2B = R_A2B.row_join(P_A2B)
    T_A2B = T_A2B.col_join(sp.Matrix([
        [0, 0, 0, 1]
    ]))

    return T_A2B

def rotation_matrix_to_rpy(rotation_matrix):
    r11 = rotation_matrix.row(0)[0]
    r21 = rotation_matrix.row(1)[0]
    r31 = rotation_matrix.row(2)[0]
    r32 = rotation_matrix.row(2)[1]
    r33 = rotation_matrix.row(2)[2]

    roll = math.atan2(r32, r33)
    pitch = math.atan2(r21, r11)
    yaw = math.atan2(-r31, math.sqrt(r11**2 + r21**2))

    sp.pprint(roll)
    sp.pprint(pitch)
    sp.pprint(yaw)

m = sp.Matrix([
    [sp.symbols("r11"), sp.symbols("r12"), sp.symbols("r13"), sp.symbols("px")],
    [sp.symbols("r21"), sp.symbols("r22"), sp.symbols("r23"), sp.symbols("py")],
    [sp.symbols("r31"), sp.symbols("r32"), sp.symbols("r33"), sp.symbols("pz")],
    [0,                 0,                 0,                 1],
])

link1 = 0.05
link2 = 0.10048
link3 = 0.1
link4 = 0.055

if __name__ == "__main__":
    list = calculate([
        DH_Param(0,             0,      0,      theta1),
        DH_Param(sp.rad(-90),   0,      d2,     theta2),
        DH_Param(0,             a2,     0,      theta3),
        DH_Param(0,             a3,     0,      theta4),
        DH_Param(sp.rad(90),    0,      d5,     theta5),
    ])

    # for T in list:
        # sp.pprint(T)

    T = product(list)
    # sp.pprint(T)

    T = T.evalf(chop=True, subs={
            a2: link2,
            a3: link3,
            theta1: math.radians(10),
            theta2: math.radians(20),
            theta3: math.radians(30),
            theta4: math.radians(40),
            theta5: math.radians(50),
            d2: link1,
            d5: link4,
            # theta1: math.radians(pulse_to_deg(500)),
            # theta2: math.radians(pulse_to_deg(600)),
            # theta3: math.radians(pulse_to_deg(820)),
            # theta4: math.radians(pulse_to_deg(110)),
            # theta5: math.radians(pulse_to_deg(500)),
        }
    )
    sp.pprint(T)
        
    r13 = T.row(0)[2]
    r23 = T.row(1)[2]
    r31 = T.row(2)[0]
    r32 = T.row(2)[1]
    r33 = T.row(2)[2]
    px = T.row(0)[3]
    py = T.row(1)[3]
    pz = T.row(2)[3]

    # Theta1
    theta_1 = math.degrees(math.atan2(r23, r13))

    # Theta5
    theta_5 = math.degrees(math.atan2(r32, -r31))

    # Theta2 + Theta3 + Theta4
    c234 = r33
    s234 = r32 / math.sin(math.radians(theta_5))
    print(c234, s234)
    theta_234 = math.degrees(math.atan2(s234, c234))

    # Theta3
    A = px*math.cos(math.radians(theta_1)) + py*math.sin(math.radians(theta_1)) - link4*math.sin(math.radians(theta_234)) # a2*c2 + a3c23 = px*c1 - py*s1 - d5*s234
    B = link4*math.cos(math.radians(theta_234)) - pz                                                                      # a2*s2 + a3s23 = d5*c234 - pz

    c3 = (A**2 + B**2 - link2**2 - link3**2)/(2*link2*link3)

    theta_3 = math.degrees(math.acos(c3)) 

    # Theta2
    C = link2 + link3*math.cos(math.radians(theta_3)) 
    D = link3*math.sin(math.radians(theta_3))
    theta_2 = math.degrees(math.atan2(B, A) - math.atan2(D, C))

    # Theta4
    theta_4 = theta_234 - theta_2 - theta_3
                           
    # print(theta_1)
    # print(theta_2)
    # print(theta_3)
    # print(theta_4)
    # print(theta_5)