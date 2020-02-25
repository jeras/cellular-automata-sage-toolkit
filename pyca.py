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
    sage: for _ in range(7) : lt.next(); lt
    '00110111110001'
    '01111100010011'
    '11000100110111'
    '01001101111100'
    '11011111000100'
    '11110001001101'
    '00010011011111'

    Creation of an image of a 1D CA run. The initial configuration is 1024
    cells with random value, the preciously defined rule 110 CA is used.
    The output image is stored into a local file "rule110_rqndom.png"

    sage: import numpy, ca_vizual
    sage: lt = pyca.lattice (ca, numpy.random.random_integers(0,1,(1024)))
    sage: ca_img = ca_vizual.array2image (lt.run(2047), 2, 'rule110_random')
"""

#import sage
import numpy
import pyca_1d     as ca1d
import pyca_common as common
#import ca_vizual

ca_neighborhood_by_name = {
    'elementary'  : ((-1,), (0,), (1,)),
    'trid'        : ((0,0), (0,1), (1,0)),
    'quad'        : ((0,0), (0,1), (1,0), (1,1)),
    'von Neumann' : ((0,-1), (-1,0), (0,0), (+1,0), (0,+1)),
    'Moore'       : ((-1,-1), (0,-1), (+1,-1), (-1,0), (0,0), (+1,0), (-1,+1), (0,+1), (+1,+1)),
    'hex'         : ()
}

class rule (object) :
    """
    A CA rule is a tuple defined by three parameters (k, m, r).
    
    INPUTS:
        k  -- number of states per cell
        a  -- neighborhood (a list of cell coordinates)
        r  -- rule number

    INTERNALS:
        m  -- neighborhood size (number of cells composing a neighborhood)
        d  -- number of dimensions
        f  -- local tansition function (lookup table)
        D  -- a list of $k$ de Bruijn matrices
        Sf -- forward state machine
        Sb -- backward state machine
    """

    def __init__(ca, k, a, r):
        ca.k = k
        ca.a = a
        ca.r = r

    def __repr__(ca):
        return "Cellular Automaton (states = "+str(ca.k)+", neighborhood = "+str(ca.a)+", rule = "+str(ca.r)+")"

    def get_a (ca) :
        return ca.__a
    def set_a (ca, a) :
        # the neighborhood is provided by its common name
        if (type(a) is str) :
            try :
                ca.__a = ca_neighborhood_by_name [a]
            except KeyError :
                print("Error: unsupported neighborhood name")
                return
        # check if the neighborhood is a tuple, otherwise translate it
        else :
            if ( (type(a   ) not in (tuple,)) or
                 (type(a[0]) not in (tuple,)) ) :
                if (type(a[0]) in (list,)) :
                    ca.__a = tuple([tuple( x ) for x in a])
                else :
                    ca.__a = tuple([tuple([x]) for x in a])
        # set neighborhood related parameters
        ca.m = len(ca.__a)
        ca.d = len(ca.__a[0])
        ca.sh = -ca.__a[0][0]
    a = property(get_a, set_a)

    def get_r (ca) :
        return ca.__r
    def set_r (ca, r) :
        ca.__r = r
        # build transition lookup table
        ca.f = numpy.array([ca.r // ca.k**n % ca.k for n in range(ca.k**ca.m)], dtype='uint8')
        # for 1D CA some rule properties can be computed
        if (ca.d == 1) :
            ca.D = [numpy.mat (numpy.zeros ((ca.k**(ca.m-1), ca.k**(ca.m-1)), int)) for k in range(ca.k)]
            for n in range(ca.k**ca.m) :
                o_l = n // ca.k; o_r = n % (ca.k**(ca.m-1))
                ca.D [ca.f[n]] [o_l, o_r] = 1
        # construct the subset diagram
        ca.Sf = [ [ common.list2int (common.list2bool (numpy.array(numpy.mat(common.int2list(i, ca.k, ca.k**(ca.m-1))) * ca.D[c])[0]), ca.k) for i in range(2**(ca.k**(ca.m-1))) ] for c in range(ca.k) ]
       #ca.Sb = [ [ common.list2int (common.list2bool (numpy.array(ca.D[c] * numpy.mat(common.int2list(i, ca.k, ca.k**(ca.m-1))))[0]), ca.k) for i in range(2**(ca.k**(ca.m-1))) ] for c in range(ca.k) ]
    r = property(get_r, set_r)

    def GoE_count(ca, N) :
        M = numpy.zeros ((2**(ca.k**(ca.m-1)), 2**(ca.k**(ca.m-1))), int)
        for c in range(ca.k) :
            for i in range(2**(ca.k**(ca.m-1))) :
                 M[i][ca.Sf[c][i]] = M[i][ca.Sf[c][i]] + 1
        M = numpy.matrix (M)
        V = numpy.zeros (2**(ca.k**(ca.m-1)), int)
        V[2**(ca.k**(ca.m-1))-1] = 1
        V = numpy.matrix (V)
        L = []
        for _ in range(N) :
            L.append (V[0,0])
            V = V*M
        #Lf = [float(L[i]/2**i) for i in range(len(L))]  # float GoE density
        return L


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
            elif (type(C) in (str,)) :
                lt.Cc = numpy.fromstring(C, dtype=dtype)-int(48)
            else :
                print("Error: incompattible configuration")
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
                    print("Error: incompattible configuration")
                    return
            else :
                print("Error: incompattible configuration")
                return
        else :
            # nD CA where n>2 can only be initialized from numpy.ndarray
            if (type(C) in (numpy.ndarray,)) :
                lt.Cc = numpy.array(C, dtype=dtype)
            else :
                print("Error: incompattible configuration")
                return
        # an empty array for neighborhoods, can be used for temporal results
        lt.Cn = numpy.empty(lt.N, dtype=dtype)
        # parsing the input boundary
        if (b in ('cyclic', 'open')) :
            lt.b = b
        else :
            print("Error: currently only 'cyclic' and 'open' boundaries are supported")
            return

    def __repr__(lt):
        return repr((lt.Cc+int(48)).tostring())

    def next (lt, t=1) :
        """
        Performs a step forward in time, the result is stored back into the
        configuration 'Cc'.
        """
        if (lt.type == '1D') :
            if (lt.b == 'cyclic') :
                for _ in range(t) :
                    ca1d.next (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc)
                shift = t * (-lt.ca.a[0][0])
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-shift:lt.N], lt.Cc[0:lt.N-shift]), axis=0)
            else :
                for _ in range(t) :
                    ca1d.next (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc)
                    lt.N = lt.N - (lt.ca.m-1)
                    lt.Cc = lt.Cc[0:lt.N]
        elif ( (lt.type == 'general') & (lt.ca.d == 2)) :
            if (lt.b == 'cyclic') :
                for _ in range(t) :
                    ca2d.ca2d_next_generic (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc)
                shift = t * (-lt.ca.a[0][0])
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-shift:lt.N], lt.Cc[0:lt.N-shift]), axis=0)
            else :
                for _ in range(t) :
                    ca2d.ca2d_next_generic (lt.ca.k, lt.ca.m, lt.ca.f, lt.N, lt.Cc)
                    lt.N = lt.N - (lt.ca.m-1)
                    lt.Cc = lt.Cc[0:lt.N]
            
        else :
            print("Error: feature not yet implemented")
        return

    def run (lt, t) :
        if (lt.ca.d == 1) :
           #if ( (lt.Cc.dtype in (dtype('uint8'),)) and (lt.ca.k == 2) ) :
            if (lt.ca.k == 2) :
                Ct = numpy.empty([t+1, lt.N], dtype='uint8')
                Ct [0] = lt.Cc
                for y in range(t) :
                    lt.next()
                    Ct [y+1] = lt.Cc
        else :
            print("Error: feature not yet implemented, only 1D CA are supported")
        return Ct

    def isGoE (lt) :
        if ((lt.ca.d == 1) and (lt.b == 'open')) :
            s = 2**(lt.ca.k**(lt.ca.m-1))-1
            for x in range(lt.N) : s = lt.ca.Sf[lt.Cc[x]][s];
            return (s == 0)
        else :
            return("Unsupported boundary")

    def prev (lt) :
        if (lt.ca.d == 1) :
         
            C_p = []

            if (lt.b == 'cyclic') :
                lt.D_x_b = [numpy.identity(lt.ca.k**(lt.ca.m-1), dtype="int")]
                for x in range(lt.N) :
                    lt.D_x_b.append (lt.ca.D[lt.Cc[lt.N-1-x]] * lt.D_x_b[x])
                lt.p = sum ([lt.D_x_b [lt.N] [i,i] for i in range(lt.ca.k**(lt.ca.m-1))])
 
                C_p = [lattice(lt.ca, lt.N*[0], lt.b) for i in range(lt.p)]
                o_p0 = [];
                for o in range(lt.ca.k**(lt.ca.m-1)) :
                    o_p0.extend([o for d in range(lt.D_x_b[lt.N][o,o])])
                o_p = list(o_p0)
               
                for x in range(lt.N) :
                    i = 0
                    while (i<lt.p) :
                        o_L = o_p[i]; o_R = o_p0[i]
                        for c in range(lt.ca.k) :
                            n = o_L * lt.ca.k + c
                            if (lt.Cc[x] == lt.ca.f[n]) :
                                o_x = n % (lt.ca.k**(lt.ca.m-1))
                                p_i = lt.D_x_b[lt.N-x-1][o_x,o_R]
                                for p_c in range(p_i) :
                                    C_p[i].Cc [(x+lt.ca.sh) % lt.N] = c
                                    o_p[i] = o_x
                                    i = i+1

                return C_p

            elif (lt.b == 'open') :
                b_L = b_R = numpy.matrix((lt.ca.k**(lt.ca.m-1))*[1])
                lt.b_x_b = [b_R.T]
#              else :
#                b_L = (lt.b[0]); b_R = vector(lt.b[1])
#                lt.b_x_b = [b_R]

                x = 0
                for x in range(lt.N) :
                    lt.b_x_b.append (lt.ca.D[lt.Cc[lt.N-1-x]] * lt.b_x_b[x])
                lt.p = (b_L * lt.b_x_b[lt.N-1])[0,0]

                C_p = [lattice(lt.ca, (lt.N+lt.ca.m-1)*[0], lt.b) for i in range(lt.p)]
                o_p = [];
                for o in range(lt.ca.k**(lt.ca.m-1)) :
                    o_p.extend([o for d in range(b_L[0,o] * lt.b_x_b[lt.N][o,0])])
                for i in range(lt.p) :
                    C_p[i].Cc [0:lt.ca.m-1] = common.int2list(o_p[i], lt.ca.k, lt.ca.m-1)

                for x in range(lt.N) :
                    i = 0
                    while (i<lt.p) :
                        o_L = o_p[i];
                        for c in range(lt.ca.k) :
                            n = o_L * lt.ca.k + c
                            if (lt.Cc[x] == lt.ca.f[n]) :
                                o_x = n % (lt.ca.k**(lt.ca.m-1))
                                p_i = lt.b_x_b[lt.N-x-1][o_x]
                                for p_c in range(p_i) :
                                    C_p[i].Cc [x+lt.ca.m-1] = c
                                    o_p[i] = o_x
                                    i = i+1

                return C_p

            else :
                print("Error: there is no preimage implementation for automata with more than 1 dimensions.")