# This function will calculate five necessary angles: theta_one, theta_two, phi_one, phi_two and omega
# Input variables: four coordinates indicating COM and interaction site of two chains

def angles(COM1, COM2, int_site1, int_site2,normal_point1,normal_point2):
    import math
    import numpy as np

    def mag(x):
        return math.sqrt(sum(i ** 2 for i in x))

    def unit(x):
        x_unit = [x[0]/mag(x), x[1]/mag(x), x[2]/mag(x)]
        return x_unit

    # calculate theta
    v1 = [int_site1[0] - COM1[0], int_site1[1] - COM1[1], int_site1[2] - COM1[2]]
    v2 = [int_site2[0] - COM2[0], int_site2[1] - COM2[1], int_site2[2] - COM2[2]]
    sigma1 = [int_site1[0] - int_site2[0], int_site1[1] - int_site2[1], int_site1[2] - int_site2[2]]
    sigma2 = [int_site2[0] - int_site1[0], int_site2[1] - int_site1[1], int_site2[2] - int_site1[2]]
    sigma_magnitude = math.sqrt(sigma1[0] ** 2 + sigma1[1] ** 2 + sigma1[2] ** 2)
    theta1 = math.acos(np.dot(v1, sigma1) / mag(v1) / mag(sigma1))
    theta2 = math.acos(np.dot(v2, sigma2) / mag(v2) / mag(sigma2))

    # calculate phi
    #n1 = [0,0,1]
    #n2 = [0,0,1]
    #if np.cross(v1, n1)[0] + np.cross(v1, n1)[1] + np.cross(v1, n1)[2] == 0:
    #    n1 = [0,1,0]
    #if np.cross(v2, n2)[0] + np.cross(v2, n2)[1] + np.cross(v2, n2)[2] == 0:
    #    n2 = [0,1,0]
    #normal_point1 = [COM1[0] + n1[0], COM1[1] + n1[1], COM1[2] + n1[2]]
    #normal_point2 = [COM2[0] + n2[0], COM2[1] + n2[1], COM2[2] + n2[2]]
    n1 = unit([normal_point1[0] - COM1[0], normal_point1[1] - COM1[1], normal_point1[2] - COM1[2]])
    n2 = unit([normal_point2[0] - COM2[0], normal_point2[1] - COM2[1], normal_point2[2] - COM2[2]])
    t1_1 = unit(np.cross(v1, sigma1))
    t2_1 = unit(np.cross(v1, n1))
    t1_2 = unit(np.cross(v2, sigma2))
    t2_2 = unit(np.cross(v2, n2))
    phi1 = math.acos(np.dot(t1_1, t2_1))
    phi2 = math.acos(np.dot(t1_2, t2_2))

    # determine the sign of phi (+/-)
    v1_uni = unit(v1)
    v2_uni = unit(v2)
    n1_proj = [n1[0] - v1_uni[0] * np.dot(v1_uni,n1), n1[1] - v1_uni[1] * np.dot(v1_uni,n1), n1[2] - v1_uni[2] * np.dot(v1_uni,n1)]
    sigma1_proj = [sigma1[0] - v1_uni[0] * np.dot(v1_uni,sigma1), sigma1[1] - v1_uni[1] * np.dot(v1_uni,sigma1), sigma1[2] - v1_uni[2] * np.dot(v1_uni,sigma1)]
    # phi1_test = math.acos(np.dot(n1_proj,sigma1_proj)/mag(n1_proj)/mag(sigma1_proj))
    n2_proj = [n2[0] - v2_uni[0] * np.dot(v2_uni, n2), n2[1] - v2_uni[1] * np.dot(v2_uni, n2),n2[2] - v2_uni[2] * np.dot(v2_uni, n2)]
    sigma2_proj = [sigma2[0] - v2_uni[0] * np.dot(v2_uni, sigma2), sigma2[1] - v2_uni[1] * np.dot(v2_uni, sigma2), sigma2[2] - v2_uni[2] * np.dot(v2_uni, sigma2)]
    # phi2_test = math.acos(np.dot(n2_proj,sigma2_proj)/mag(n2_proj)/mag(sigma2_proj))
    phi1_dir = unit(np.cross(sigma1_proj, n1_proj))
    phi2_dir = unit(np.cross(sigma2_proj, n2_proj))
    if abs(v1_uni[0] - phi1_dir[0]) < 10 ** -10:
        phi1 = -phi1
    elif abs(v1_uni[0] + phi1_dir[0]) < 10 ** -10:
        phi1 = phi1
    else:
        print("Wrong phi1 angle.")
    if abs(v2_uni[0] - phi2_dir[0]) < 10 ** -10:
        phi2 = -phi2
    elif abs(v2_uni[0] + phi2_dir[0]) < 10 ** -10:
        phi2 = phi2
    else:
        print("Wrong phi2 angle.")

    # calculate omega
    a1 = np.cross(sigma1, v1) / mag(np.cross(sigma1, v1))
    a2 = np.cross(sigma1, v2) / mag(np.cross(sigma1, v2))
    omega = math.acos(np.dot(a1, a2))

    # determine the sign of omega (+/-)
    sigma1_uni = unit(sigma1)
    v1_proj = [v1[0] - sigma1_uni[0] * np.dot(sigma1_uni, v1), v1[1] - sigma1_uni[1] * np.dot(sigma1_uni, v1),
               v1[2] - sigma1_uni[2] * np.dot(sigma1_uni, v1)]
    v2_proj = [v2[0] - sigma1_uni[0] * np.dot(sigma1_uni, v2), v2[1] - sigma1_uni[1] * np.dot(sigma1_uni, v2),
               v2[2] - sigma1_uni[2] * np.dot(sigma1_uni, v2)]
    # omega_test = math.acos(np.dot(v1_proj, v2_proj) / mag(v2_proj) / mag(v1_proj))
    omega_dir = unit(np.cross(v1_proj, v2_proj))
    if abs(sigma1_uni[0] - omega_dir[0]) < 10 ** -10:
        omega = -omega
    elif abs(sigma1_uni[0] + omega_dir[0]) < 10 ** -10:
        omega = omega
    else:
        print("Wrong omega angle.")

    return theta1, theta2, phi1, phi2, omega, sigma_magnitude, normal_point1, normal_point2
