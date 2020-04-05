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
    sage: ca.transition_table
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
    A CA rule is a tuple defined by three parameters (stateset, neighborhood, rule).
    
    INPUTS:
        stateset     -- number of states per cell
        neighborhood -- a list of cell coordinates
        rule         -- rule number

    INTERNALS:
        neighborhood_size -- number of cells composing a neighborhood
        dimensions        -- number of dimensions
        transiotion_table -- local tansition function (lookup table)
        D  -- a list of $stateset$ de Bruijn matrices
        Sf -- forward state machine
        Sb -- backward state machine
    """

    def __init__(ca, stateset, neighborhood, rule):
        ca.stateset     = stateset
        ca.neighborhood = neighborhood
        ca.rule         = rule

    def __repr__(ca):
        return "Cellular Automaton (stateset = "+str(ca.stateset)+", neighborhood = "+str(ca.neighborhood)+", rule = "+str(ca.rule)+")"

    def get_neighborhood (ca) :
        return ca.__neighborhood
    def set_neighborhood (ca, neighborhood) :
        # the neighborhood is provided by its common name
        if (type(neighborhood) is str) :
            try :
                ca.__neighborhood = ca_neighborhood_by_name [neighborhood]
            except KeyError :
                print("Error: unsupported neighborhood name")
                return
        # check if the neighborhood is a tuple, otherwise translate it
        else :
            if ( (type(neighborhood   ) not in (tuple,)) or
                 (type(neighborhood[0]) not in (tuple,)) ) :
                if (type(neighborhood[0]) in (list,)) :
                    ca.__neighborhood = tuple([tuple( x ) for x in neighborhood])
                else :
                    ca.__neighborhood = tuple([tuple([x]) for x in neighborhood])
        # set neighborhood related parameters
        ca.neighborhood_size = len(ca.__neighborhood)
        ca.dimensions        = len(ca.__neighborhood[0])
        ca.sh = -ca.__neighborhood[0][0]
    neighborhood = property(get_neighborhood, set_neighborhood)

    def get_rule (ca) :
        return ca.__rule
    def set_rule(ca, rule) :
        ca.__rule = rule
        # build transition lookup table
        ca.transition_table = numpy.array([ca.rule // ca.stateset**n % ca.stateset for n in range(ca.stateset**ca.neighborhood_size)], dtype='uint8')
        # for 1D CA some rule properties can be computed
        if (ca.dimensions == 1) :
            ca.de_Bruijn_matrices = [numpy.mat (numpy.zeros ((ca.stateset**(ca.neighborhood_size-1), ca.stateset**(ca.neighborhood_size-1)), int)) for k in range(ca.stateset)]
            for n in range(ca.stateset**ca.neighborhood_size) :
                o_l = n // ca.stateset; o_r = n % (ca.stateset**(ca.neighborhood_size-1))
                ca.de_Bruijn_matrices [ca.transition_table[n]] [o_l, o_r] = 1
        # construct the subset diagram
        ca.Sf = [ [ common.list2int (common.list2bool (numpy.array(numpy.mat(common.int2list(i, ca.stateset, ca.stateset**(ca.neighborhood_size-1))) * ca.de_Bruijn_matrices[c])[0]), ca.stateset) for i in range(2**(ca.stateset**(ca.neighborhood_size-1))) ] for c in range(ca.stateset) ]
       #ca.Sb = [ [ common.list2int (common.list2bool (numpy.array(ca.de_Bruijn_matrices[c] * numpy.mat(common.int2list(i, ca.stateset, ca.stateset**(ca.neighborhood_size-1))))[0]), ca.stateset) for i in range(2**(ca.stateset**(ca.neighborhood_size-1))) ] for c in range(ca.stateset) ]
    rule = property(get_rule, set_rule)

    def GoE_count(ca, N) :
        M = numpy.zeros ((2**(ca.stateset**(ca.neighborhood_size-1)), 2**(ca.stateset**(ca.neighborhood_size-1))), int)
        for c in range(ca.stateset) :
            for i in range(2**(ca.stateset**(ca.neighborhood_size-1))) :
                 M[i][ca.Sf[c][i]] = M[i][ca.Sf[c][i]] + 1
        M = numpy.matrix (M)
        V = numpy.zeros (2**(ca.stateset**(ca.neighborhood_size-1)), int)
        V[2**(ca.stateset**(ca.neighborhood_size-1))-1] = 1
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
        configuration -- CA configuration (can be a list, n-D array or string)
        boundary      -- boundary type (default 'cyclic')
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
    def __init__(lt, ca, configuration, boundary='cyclic', dtype='uint8') :
        lt.ca = ca
        lt.N  = len(configuration)
        # detecting the CA type
        if ( (lt.ca.dimensions==1) & (lt.ca.neighborhood == tuple([tuple([x+lt.ca.neighborhood[0][0]]) for x in range(lt.ca.neighborhood_size)])) ) :
            lt.type = '1D'
        else :
            lt.type = 'general'
        # parsing the input configuration
        if (lt.ca.dimensions == 1) :
            # 1D CA can be initialized from numpy.ndarray, list or string
            if (type(configuration) in (numpy.ndarray, list)) :
                lt.Cc = numpy.array(configuration, dtype=dtype)
            elif (type(configuration) in (str,)) :
                lt.Cc = numpy.fromstring(configuration, dtype=dtype)-int(48)
            else :
                print("Error: incompattible 1D configuration")
                return
        elif (lt.ca.dimensions == 2) :
            # 2D CA can be initialized from numpy.ndarray, list of lists or strings, or from images
            if (type(configuration) in (numpy.ndarray,)) :
                lt.Cc = numpy.array(C, dtype=dtype)
            elif (type(configuration) in (list,)) :
                if (type(configuration[0]) in (list,)) :
                    lt.Cc = numpy.array(configuration, dtype=dtype)
                elif (type(configuration) in (string,)) :
                    lt.Cc = numpy.empty(lt.N, dtype=dtype)
                    for y in range(len(configuration)) :
                        lt.Cc = numpy.fromstring(configuration[y], dtype=dtype)-int(48)
                else :
                    print("Error: incompattible 2D configuration")
                    return
            else :
                print("Error: incompattible 2D configuration")
                return
        else :
            # nD CA where n>2 can only be initialized from numpy.ndarray
            if (type(configuration) in (numpy.ndarray,)) :
                lt.Cc = numpy.array(configuration, dtype=dtype)
            else :
                print("Error: incompattible nD configuration")
                return
        # an empty array for neighborhoods, can be used for temporal results
        lt.Cn = numpy.empty(lt.N, dtype=dtype)
        # parsing the input boundary
        if (boundary in ('cyclic', 'open')) :
            lt.boundaryoundary = boundary
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
            if (lt.boundary == 'cyclic') :
                for _ in range(t) :
                    ca1d.next (lt.ca.stateset, lt.ca.neighborhood_size, lt.ca.transition_table, lt.N, lt.Cc)
                shift = t * (-lt.ca.neighborhood[0][0])
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-shift:lt.N], lt.Cc[0:lt.N-shift]), axis=0)
            else :
                for _ in range(t) :
                    ca1d.next (lt.ca.stateset, lt.ca.neighborhood_size, lt.ca.transition_table, lt.N, lt.Cc)
                    lt.N = lt.N - (lt.ca.neighborhood_size-1)
                    lt.Cc = lt.Cc[0:lt.N]
        elif ( (lt.type == 'general') & (lt.ca.dimensions == 2)) :
            if (lt.boundary == 'cyclic') :
                for _ in range(t) :
                    ca2d.ca2d_next_generic (lt.ca.stateset, lt.ca.neighborhood_size, lt.ca.transition_table, lt.N, lt.Cc)
                shift = t * (-lt.ca.neighborhood[0][0])
                lt.Cc = numpy.concatenate((lt.Cc[lt.N-shift:lt.N], lt.Cc[0:lt.N-shift]), axis=0)
            else :
                for _ in range(t) :
                    ca2d.ca2d_next_generic (lt.ca.stateset, lt.ca.neighborhood_size, lt.ca.transition_table, lt.N, lt.Cc)
                    lt.N = lt.N - (lt.ca.neighborhood_size-1)
                    lt.Cc = lt.Cc[0:lt.N]
            
        else :
            print("Error: feature not yet implemented")
        return

    def run (lt, t) :
        if (lt.ca.dimensions == 1) :
           #if ( (lt.Cc.dtype in (dtype('uint8'),)) and (lt.ca.stateset == 2) ) :
            if (lt.ca.stateset == 2) :
                Ct = numpy.empty([t+1, lt.N], dtype='uint8')
                Ct [0] = lt.Cc
                for y in range(t) :
                    lt.next()
                    Ct [y+1] = lt.Cc
        else :
            print("Error: feature not yet implemented, only 1D CA are supported")
        return Ct

    def isGoE (lt) :
        if ((lt.ca.dimensions == 1) and (lt.boundary == 'open')) :
            s = 2**(lt.ca.stateset**(lt.ca.neighborhood_size-1))-1
            for x in range(lt.N) : s = lt.ca.Sf[lt.Cc[x]][s];
            return (s == 0)
        else :
            return("Unsupported boundary")

    def prev (lt) :
        if (lt.ca.dimensions == 1) :
         
            C_p = []

            if (lt.boundary == 'cyclic') :
                lt.D_x_b = [numpy.identity(lt.ca.stateset**(lt.ca.neighborhood_size-1), dtype="int")]
                for x in range(lt.N) :
                    lt.D_x_b.append (lt.ca.de_Bruijn_matrices[lt.Cc[lt.N-1-x]] * lt.D_x_b[x])
                lt.p = sum ([lt.D_x_b [lt.N] [i,i] for i in range(lt.ca.stateset**(lt.ca.neighborhood_size-1))])
 
                C_p = [lattice(lt.ca, lt.N*[0], lt.boundary) for i in range(lt.p)]
                o_p0 = [];
                for o in range(lt.ca.stateset**(lt.ca.neighborhood_size-1)) :
                    o_p0.extend([o for d in range(lt.D_x_b[lt.N][o,o])])
                o_p = list(o_p0)
               
                for x in range(lt.N) :
                    i = 0
                    while (i<lt.p) :
                        o_L = o_p[i]; o_R = o_p0[i]
                        for c in range(lt.ca.stateset) :
                            n = o_L * lt.ca.stateset + c
                            if (lt.Cc[x] == lt.ca.transition_table[n]) :
                                o_x = n % (lt.ca.stateset**(lt.ca.neighborhood_size-1))
                                p_i = lt.D_x_b[lt.N-x-1][o_x,o_R]
                                for p_c in range(p_i) :
                                    C_p[i].Cc [(x+lt.ca.sh) % lt.N] = c
                                    o_p[i] = o_x
                                    i = i+1

                return C_p

            elif (lt.boundary == 'open') :
                b_L = b_R = numpy.matrix((lt.ca.stateset**(lt.ca.neighborhood_size-1))*[1])
                lt.boundary_x_b = [b_R.T]
#              else :
#                b_L = (lt.boundary[0]); b_R = vector(lt.boundary[1])
#                lt.boundary_x_b = [b_R]

                x = 0
                for x in range(lt.N) :
                    lt.boundary_x_b.append (lt.ca.de_Bruijn_matrices[lt.Cc[lt.N-1-x]] * lt.boundary_x_b[x])
                lt.p = (b_L * lt.boundary_x_b[lt.N-1])[0,0]

                C_p = [lattice(lt.ca, (lt.N+lt.ca.neighborhood_size-1)*[0], lt.boundary) for i in range(lt.p)]
                o_p = [];
                for o in range(lt.ca.stateset**(lt.ca.neighborhood_size-1)) :
                    o_p.extend([o for d in range(b_L[0,o] * lt.boundary_x_b[lt.N][o,0])])
                for i in range(lt.p) :
                    C_p[i].Cc [0:lt.ca.neighborhood_size-1] = common.int2list(o_p[i], lt.ca.stateset, lt.ca.neighborhood_size-1)

                for x in range(lt.N) :
                    i = 0
                    while (i<lt.p) :
                        o_L = o_p[i];
                        for c in range(lt.ca.stateset) :
                            n = o_L * lt.ca.stateset + c
                            if (lt.Cc[x] == lt.ca.transition_table[n]) :
                                o_x = n % (lt.ca.stateset**(lt.ca.neighborhood_size-1))
                                p_i = lt.boundary_x_b[lt.N-x-1][o_x]
                                for p_c in range(p_i) :
                                    C_p[i].Cc [x+lt.ca.neighborhood_size-1] = c
                                    o_p[i] = o_x
                                    i = i+1

                return C_p

            else :
                print("Error: there is no preimage implementation for automata with more than 1 dimensions.")