import numpy as np

def make_quat_rot1( x3):
    import math
    pi = math.acos(-1.0)
    x0 = x3[0];x1 = x3[1] ; x2 = x3[2]
    quat_rot = [(1.0 - x0)**0.5 * math.sin(2.0*pi*x1) ,(1.0 - x0)**0.5 * math.cos(2.0*pi*x1) ,x0**0.5 * math.sin(2.0*pi*x2) , x0**0.5 * math.cos(2.0*pi*x2)]
    return np.array(quat_rot)

def quaternion_mult ( q1, q2):
    q3 = np.zeros(4)#range(0,4)
    q3[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    q3[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]
    q3[2] = q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3]
    q3[3] = q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1]
    return np.array(q3)

def quaternion_rotatn(vec,quat):
#    dhc- this has a silly name,, but the old routine with a similar name is still used by pete
    quat_inv = np.array([quat[0],-1.0*quat[1],-1.0*quat[2],-1.0*quat[3]])
    return quaternion_mult( quaternion_mult(quat, np.array([0.0, vec[0], vec[1], vec[2]])), quat_inv )[1:4]
