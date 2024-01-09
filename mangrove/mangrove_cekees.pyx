#cython: language_level=3
import numpy as np
#from libc.math cimport csqrtf
#from libc.math cimport cos, sin, sqrt, fmin, fabs, csqrt, cabs
from math import cos, sin, sqrt, sqrt
from cmath import sqrt as csqrt
cimport numpy as np
cimport cython
np.import_array()


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def sdf_vectorized(double t, np.ndarray x,
                   np.ndarray phi):
    cdef np.ndarray [np.float64_t, ndim=1] L = np.array([3.22,1.0,1.0])
    cdef np.ndarray[np.float64_t, ndim=1] theta=np.array(
        [0.0, 0.7853981633974483, 1.5707963267948966, 2.356194490192345, 3.141592653589793, 3.9269908169872414, 4.71238898038469]) 
    cdef np.ndarray[np.float64_t, ndim=1] X_r=np.array(
        [0.1595238095238095, 0.21666666666666667, 0.2714285714285714, 0.3285714285714285, 0.38571428571428573, 0.44285714285714284, 0.5])
    cdef np.ndarray[np.float64_t, ndim=1] H_r=np.array(
        [0.13809523809523808, 0.16904761904761903, 0.19999999999999998, 0.23095238095238094, 0.25952380952380955, 0.29047619047619044, 0.32142857142857145])
    cdef int n = phi.shape[0]
    cdef int nr = theta.shape[0]
    cdef int ti, i, ri, gri
    cdef double R,X_0,Y_0,Z_0,A,B,dir_x, dir_y
    cdef double G
    cdef np.ndarray[np.float64_t, ndim=1] tx = np.array([2.76,3.04,3.32,3.6,3.88,2.76,3.04,3.32,3.6,3.88,2.62,2.9,3.18,3.46,3.74,4.02]) 
    cdef np.ndarray[np.float64_t, ndim=1] ty = np.array([0.0625,0.0625,0.0625,0.0625,0.0625,0.2375,0.2375,0.2375,0.2375,0.2375,0.15,0.15,0.15,0.15,0.15,0.15])
    cdef np.ndarray[np.complex128_t, ndim=1] G_root = np.zeros((3,), dtype=np.complex128)
    cdef double complex num, denom
    for i in range(n):
        phi[i] = 1.0e16
        for ti in range(len(tx)):
            X_0 = x[i,0] - tx[ti]
            Y_0 = x[i,1] - ty[ti]
            Z_0 = x[i,2]
            phi[i] = min(phi[i], sqrt(X_0**2 + Y_0**2) - 0.044)  #0.02721428571428571)
#            for ri in range(nr):
#                dir_x  = cos(theta[ri]+(ti+1)*2.0943951023931953)#perturb by 2pi/3 so the middle tree has the referene orientation
#                dir_y  = sin(theta[ri]+(ti+1)*2.0943951023931953)
#                A = -H_r[ri]/X_r[ri]**2
#                B = H_r[ri]
#                num = 2**(1/3)*3**(2/3)*(3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3)*(1 + csqrt(3)*1.0j)**2 + 12*2**(1/3)*(-2*A*B + 2*A*Z_0 - dir_x**2 - dir_y**2))
#                denom = 36*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)*(1 + csqrt(3)*1.0j)
#                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
#                    G_root[0] = num/denom
#                else:
#                    G_root[0] = num
#                num = 2**(1/3)*3**(2/3)*(3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3)*(1 - csqrt(3)*1.0j)**2 + 12*2**(1/3)*(-2*A*B + 2*A*Z_0 - dir_x**2 - dir_y**2))
#                denom = 36*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)*(1 - csqrt(3)*1.0j)
#                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
#                    G_root[1] = num/denom
#                else:
#                    G_root[1] = num
#                num = 2**(1/3)*3**(2/3)*(-3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3) + 3*2**(1/3)*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2))
#                denom = 18*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)
#                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
#                    G_root[2] = num/denom
#                else:
#                    G_root[2] = num
#                for gri in range(3):
#                    R = G_root[gri].real
#                    if R < -X_r[ri]:
#                        R = -X_r[ri]
#                    if R > X_r[ri]:
#                        R = X_r[ri]
#                    G = sqrt((R*dir_x - X_0)**2 + (R*dir_y - Y_0)**2 + (A*R**2 + B - Z_0)**2) - 0.00680952380952381;
#                    phi[i] = min(phi[i], G)
