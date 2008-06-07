"""
SAGE/Python Cellular Automata Toolkit Module

Toolkit for theretical and prectical research of cellular automata. The
theoretical part focuses on the computation of preimages, the practical on
creating images and videos of CA. The proposed pronounciation of the "pyca"
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
        b  -- boundary type (the default 'cyclic' or a configuration of length ca.m-1)
        dtype -- cell value data type (default is 'uint8')

    INTERNALS:
        N  -- lattice size as number of cels in each dimension
        Cc -- internal 'uint8' array representation of the cell lattice
        Cn -- internal 'uint8' array representation of the neighborhood lattice

    Since internal array tipes are of size 'uint8', the number of states per cell
    is limited to 256.

    If the CA boundary is defined as 'cyclic' than actually there is no boundary,
    otherwise the boundary can be defined explicitely as a concatenation of the
    left and the right boundary [b_L+b_R], which can be of the same types as the
    configuration.
    """
    def __init__(lt, ca, C, b='cyclic', dtype='uint8') :
        lt.ca = ca
        lt.N  = len(C)
        # detecting the CA type
        if ( (lt.ca.d==1) & (lt.ca.a == tuple([tuple([x+lt.ca.a[0][0]]) for x in range(lt.ca.m)])) ) :
            lt.type = '1D'
            lt.sh = - lt.ca.a[0][0]
        else :
            lt.type = 'general'
        # parsing the input configuration
        if (type(C) in (numpy.ndarray, list)) :
            lt.Cc = numpy.array(C, dtype=dtype)
        else :
            dx = [0 for _ in range(lt.ca.d)]
            # TODO generalize for any number of dimensions
            lt.Cc = numpy.fromstring(C, dtype=dtype)-int(48)
        lt.Cn = numpy.empty(lt.N, dtype=dtype)
        # parsing the input boundary
        if (b == 'cyclic') :
            lt.b = b
        else :
            if (type(b) in (numpy.ndarray, list)) :
                lt.b = numpy.array(b, dtype=dtype)
            else :
                lt.b = numpy.fromstring(b, dtype=dtype)-int(48)

    def __repr__(lt):
        return repr((lt.Cc+int(48)).tostring())

    def next (lt) :
        """
        Performs a step forward in time, the result is stored back into the
        configuration 'Cc'. As a partial rasult the neighborhood configuration is
        computed an stored into 'Cn'.
        """
        if (lt.type == '1D') :
            ca1d.ca1d_next_generic (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc, lt.Cn)
            if (lt.b == 'cyclic') :
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-lt.sh:lt.N], lt.Cc[0:lt.N-lt.sh]), axis=0)
            else :
                lt.Cc = lt.Cc[0:lt.N-(ca.m-1)]
            return

