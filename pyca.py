"""
SAGE/Python Cellular Automata Toolkit Module

Toolkit for theoretical and practical research of cellular automata. The
theoretical part focuses on the computation of preimages, the practical on
creating images and videos of CA. The proposed pronunciation of the "pyca"
toolkit is the same as the Italian dish pizza.

AUTHORS:
    -- Iztok Jeras

EXAMPLES:
    
    Creation of a CA rule object:

    sage: import pyca
    sage: ca = pyca.rule (2, (-1,0,1), 110)
    sage: ca
    Cellular Automaton (states = 2, neighborhood = ((-1,), (0,), (1,)), rule = 110)
    sage: ca.f
    array([0, 1, 1, 1, 0, 1, 1, 0], dtype=uint8)

    Creation of 1D lattice, three different data types are used for the
    configuration, the results are checked to be equal:

    sage: lt = pyca.lattice (ca, [0,0,0,1,0,0,1,1,0,1,1,1,1,1])
    sage: lt
    '00010011011111'
    sage: for _ in xrange(7) : lt.next(); lt
    ....:
    '00110111110001'
    '01111100010011'
    '11000100110111'
    '01001101111100'
    '11011111000100'
    '11110001001101'
    '00010011011111'
"""

import numpy
import Image
import ca1d


class rule () :
    """
    A CA rule is a tuple defined by three parameters (k, m, r).
    
    INPUTS:
        k  -- number of states per cell
        a  -- neighborhood (a list of cell coordinates)
        r  -- rule number

    INTERNALS:
        m  -- neighborhood size (number of cells composing a neighborhood)
        d  -- number of dimensions
        D  -- a list of $k$ de Bruijn matrices
        Sf -- forward state machine
        Sb -- backward state machine
    """
    def __init__(ca, k, a, r):
        ca.k = k
        ca.a = a
        ca.m = len(ca.a)
       #if (type(ca.a[0]) not in (tuple, list))
        if ( (type(ca.a   ) not in (tuple,)) or
             (type(ca.a[0]) not in (tuple,)) ) :
           #ca.a = [[x] for x in ca.a]
            ca.a = tuple([tuple([x]) for x in ca.a])
        ca.d = len(ca.a[0])
        ca.r = r
        ca.f = numpy.array([ca.r // ca.k**n % ca.k for n in xrange(ca.k**ca.m)], dtype='uint8')

        # for 1D CA some rule properties can be computed
        if (ca.d == 1) :
            ca.D = [numpy.mat (numpy.zeros ((ca.k**(ca.m-1), ca.k**(ca.m-1)), int)) for k in xrange(ca.k)]
            for n in xrange(ca.k**ca.m) :
                o_l = n // ca.k; o_r = n % (ca.k**(ca.m-1))
                ca.D [ca.f[n]] [o_l, o_r] = 1

    #        ca.Sf = [ [ list2int (list2bool (vector(int2list(i, ca.k, ca.k**(ca.m-1))) * ca.D[c]), ca.k) for i in xrange(2**(ca.k**(ca.m-1))) ] for c in xrange(ca.k) ]
    #        ca.Sb = [ [ list2int (list2bool (ca.D[c] * vector(int2list(i, ca.k, ca.k**(ca.m-1)))), ca.k) for i in xrange(2**(ca.k**(ca.m-1))) ] for c in xrange(ca.k) ]

    def __repr__(ca):
        return "Cellular Automaton (states = "+str(ca.k)+", neighborhood = "+str(ca.a)+", rule = "+str(ca.r)+")"


class lattice () :
    """
    CA lattice

    INPUTS:
        ca -- CA rule object
        C  -- CA configuration (can be a list, n-D array or string)
        b  -- boundary type (default 'cyclic')
        dtype -- cell value data type (default 'uint8')

    METHODS:
        next(t) -- compute t time steps (default t is 1)

    INTERNALS:
        N  -- lattice size as number of cells in each dimension
        Cc -- internal array representation of the cell lattice
        Cn -- internal array representation of the neighborhood lattice

    The configuration can be initialized from different data types depending on
    the number of dimensions:
    1D -- numpy.ndarray, list or string
    2D -- numpy.ndarray, list of lists or strings, or from images (png, bmp, ...)
    nD -- numpy.ndarray only

    By default the internal arrays are of type 'uint8', the number of states
    per cell is limited to 256.

    The supported boundary tipes are 'cyclic and 'open'.
    The 'cyclic' or periodic boundary is used when the each boundary is
    directly conected to the oposite boundary (1D - start <=> end,
    2D - left <=> right, top <=> bottom). 
    The 'open' boundary is used when the defined configuration is only a
    segment of an infinite configuration, where cells outside the observed
    segment can have any value.

    The next(t) method computes t discrete time steps, the prev() method
    computes a list of preimages, the results depends on the boundary type.
    """
    def __init__(lt, ca, C, b='cyclic', dtype='uint8') :
        lt.ca = ca
        lt.N  = len(C)
        # detecting the CA type
        if ( (lt.ca.d==1) & (lt.ca.a == tuple([tuple([x+lt.ca.a[0][0]]) for x in range(lt.ca.m)])) ) :
            lt.type = '1D'
        else :
            lt.type = 'general'
        # parsing the input configuration
        if (lt.ca.d == 1) :
            # 1D CA can be initialized from numpy.ndarray, list or string
            if (type(C) in (numpy.ndarray, list)) :
                lt.Cc = numpy.array(C, dtype=dtype)
            elif (type(C) in (string,)) :
                lt.Cc = numpy.fromstring(C, dtype=dtype)-int(48)
            else :
                print "Error: incompattible configuration"
                return
        elif (lt.ca.d == 1) :
            # 2D CA can be initialized from numpy.ndarray, list of lists or strings, or from images
            if (type(C) in (numpy.ndarray,)) :
                lt.Cc = numpy.array(C, dtype=dtype)
            elif (type(C) in (list,)) :
                if (type(C[0]) in (list,)) :
                    lt.Cc = numpy.array(C, dtype=dtype)
                elif (type(C) in (string,)) :
                    lt.Cc = numpy.empty(lt.N, dtype=dtype)
                    for y in range(len(C)) :
                        lt.Cc = numpy.fromstring(C[y], dtype=dtype)-int(48)
                else :
                    print "Error: incompattible configuration"
                    return
            else :
                print "Error: incompattible configuration"
                return
        else :
            # nD CA where n>2 can only be initialized from numpy.ndarray
            if (type(C) in (numpy.ndarray,)) :
                lt.Cc = numpy.array(C, dtype=dtype)
            else :
                print "Error: incompattible configuration"
                return
        # an empty array for neighborhoods, can be used for temporal results
        lt.Cn = numpy.empty(lt.N, dtype=dtype)
        # parsing the input boundary
        if (b in ('cyclic', 'open')) :
            lt.b = b
        else :
            print "Error: currently only 'cyclic' and 'open' boundaries are supported"
            return

    def __repr__(lt):
        return repr((lt.Cc+int(48)).tostring())

    def next (lt, t=1) :
        """
        Performs a step forward in time, the result is stored back into the
        configuration 'Cc'. As a partial result the neighborhood configuration is
        computed an stored into 'Cn'.
        """
        if (lt.type == '1D') :
            if (lt.b == 'cyclic') :
                for _ in xrange(t) :
                    ca1d.ca1d_next_generic (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc, lt.Cn)
                shift = t * (-lt.ca.a[0][0])
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-shift:lt.N], lt.Cc[0:lt.N-shift]), axis=0)
            else :
                for _ in xrange(t) :
                    ca1d.ca1d_next_generic (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc, lt.Cn)
                    lt.N = lt.N - (lt.ca.m-1)
                    lt.Cc = lt.Cc[0:lt.N]
        else :
            print "Error: feature not yet implemented"
        return

    def cfg2image (lt, name) :
        """
        writes the current configuration into an image file
        """
        return

    def run2image (lt, t, name) :
        """
        computes t time steps and writes the resulting configurations into an image or video
        """
        if (lt.ca.d == 1) :
           #if ( (lt.Cc.dtype in (dtype('uint8'),)) and (lt.ca.k == 2) ) :
            if (lt.ca.k == 2) :
                cc_t = numpy.empty([t+1, lt.N], dtype='uint8')
                cc_t [0] = lt.Cc
                for y in range(t) :
                    lt.next()
                    cc_t [y+1] = lt.Cc
                img = Image.fromstring('P', (cc_t.shape[1], cc_t.shape[0]), cc_t.tostring());
                img.putpalette([0,0,0,255,255,255]);
                img.save(name+"%04u.png"%t, "PNG")
            else :
                print "Error: feature not yet implemented, unsuported 1D CA parameters"
        else :
            print "Error: feature not yet implemented"
        return

