"""
library for 1D CA performance

For documentation see the pyca module documentation.
"""

cdef extern from 'arrayobject.h' :
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

#cdef public ca1d_next_generic (int k, int m, ndarray f, ndarray Cc, ndarray Cn, int N, int sh, int b) :
def ca1d_next_generic (int k, int m, ndarray f, int N, ndarray Cc, ndarray Cn) :
    """
    Generic version of 1D CA processor.

    INPUTS:
        int k -- number of states per cell
        int m -- number of cells inside a neighborhood
        ndarray f  -- local transition function
        int N -- configuration size
        ndarray Cc -- cell configuration
        ndarray Cn -- neighborhood configuration (nut used as input, but overwritten)
    """

    cdef unsigned char *pf = <unsigned char *>f.data
    cdef unsigned char *pc = <unsigned char *>Cc.data
    cdef unsigned char *pn = <unsigned char *>Cn.data
    cdef int x  # current position inside the configuration
    cdef int n  # current neighborhood value
    cdef int w  # the weight of the first cell inside the neighborhood 

    if 1 :
        # initialization
        n = 0
        w = 1
        for x from 0 <= x < m-1 :
            n = n * k + pc[x]  # preparing the first overlap
            w = w * k          # seting the weight
     
        # building of neighborhoods
        for x from 0 <= x < N :
            n = n * k + pc[(x+m-1)%N]
            pn[x] = n
            n = n - pc[x] * w

    else :
        # initialization
        n = 0
        w = 1
        for x from 0 <= x < m-1 :
            n = n << 1 | pc[x]  # preparing the first overlap
            w = w << 1
        w = w-1

        # building of neighborhoods
        for x from 0 <= x < N :
            n = n << 1 | pc[x+m-1]
            pn[x] = n
            n = n & w

    # applying transition function
    for x from 0 <= x < N :
        pc[x] = pf[pn[x]]

    return

