#
# library for 2D CA performance
#

import sage.rings.integer

cdef extern from "arrayobject.h":
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

#
# 2D CA evolution
#

def ca2d_next_generic (int k, int m, ndarray a, ndarray f, ndarray Cc, ndarray Cn, int Nx, int Ny) :
    cdef unsigned char *pa = <unsigned char *>a.data
    cdef unsigned char *pf = <unsigned char *>f.data
    cdef unsigned char *pc = <unsigned char *>Cc.data
    cdef unsigned char *pn = <unsigned char *>Cn.data
    cdef unsigned int x, y, xp, yp, i, n

    # building of neighborhoods
    for y from 0 <= y < Ny :
        for x from 0 <= x < Nx :
            n = 0
            for i from 0 <= i < m :
                yp = y + pa[2*i]
                if (yp>=Ny) : yp = yp-Ny
                xp = x + pa[2*i+1]
                if (xp>=Nx) : xp = xp-Nx
                n = n * k + pc[yp*Nx + xp]
            pn[y*Nx + x] = n
            pc[y*Nx + x] = pf[n]

    return

def ca2d_next_fast (int k, int m, ndarray a, ndarray f, ndarray Cc, ndarray Cn, int Nx, int Ny) :
    cdef unsigned char *pa = <unsigned char *>a.data
    cdef unsigned char *pf = <unsigned char *>f.data
    cdef unsigned char *pc = <unsigned char *>Cc.data
    cdef unsigned char *pn = <unsigned char *>Cn.data
    cdef unsigned int x, y, xp, yp, xm, ym, xs, ys, i, n

    # creating masks
    xm = Nx-1; xs = 7
    ym = Ny-1

    # building of neighborhoods
    for y from 0 <= y <= ym :
        for x from 0 <= x <= xm :
            n = 0
            for i from 0 <= i < m :
                yp = y + pa[2*i]
                yp = yp & ym
                xp = x + pa[2*i+1]
                xp = xp & xm
                n = n * k + pc[yp<<xs + xp]
            pn[y<<xs + x] = n
            pc[y<<xs + x] = pf[n]

    return

