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
    Cellular Automaton (states = 2, neighborhood = [[-1,], [0,], [1,]], rule = 110)
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

    sage: import numpy as np
    sage: import ca_vizual
    sage: lt = pyca.lattice (ca, np.random.random_integers(0,1,(1024)))
    sage: ca_img = ca_vizual.array2image (lt.run(2047), 2, 'rule110_random')
"""

#import sage
import numpy       as np
#import ca_vizual

ca_neighborhood_by_name = {
    'elementary'  : [[-1,], [0,], [1,]],
    'trid'        : [[ 0, 0], [ 0,+1],
                     [+1, 0]         ],
    'quad'        : [[ 0, 0], [ 0,+1],
                     [+1, 0], [+1,+1]],
    'von Neumann' : [         [ 0,-1],
                     [-1, 0], [ 0, 0], [+1, 0],
                              [ 0,+1]         ],
    'Moore'       : [[-1,-1], [ 0,-1], [+1,-1],
                     [-1, 0], [ 0, 0], [+1, 0],
                     [-1,+1], [ 0,+1], [+1,+1]],
    'hex'         : []
}

class rule(object):
    """
    A CA rule is a tuple defined by three parameters (stateset, neighborhood, rule).
    
    INPUTS:
        stateset     -- number of states per cell
        neighborhood -- a list of cell coordinates
        rule         -- rule number

    INTERNALS:
        size -- number of cells composing a neighborhood
        dimensions        -- number of dimensions
        transiotion_table -- local tansition function (lookup table)
        D  -- a list of $stateset$ de Bruijn matrices
        Sf -- forward state machine
        Sb -- backward state machine
    """

    def __init__(ca, stateset, neighborhood, rule, dtype=np.uint8):
        ca.dtype        = dtype
        ca.stateset     = np.uint(stateset)
        ca.neighborhood = neighborhood
        ca.rule         = rule

    def __repr__(ca):
        return "Cellular Automaton (stateset = "+str(ca.stateset)+", neighborhood = "+str(ca.neighborhood)+", rule = "+str(ca.rule)+", dtype = "+repr(ca.dtype)+")"

    def get_neighborhood(ca):
        return ca.__neighborhood

    def set_neighborhood(ca, neighborhood):
        # the neighborhood is provided by its common name
        if (type(neighborhood) is str) :
            try :
                ca.__neighborhood = ca_neighborhood_by_name [neighborhood]
            except KeyError :
                print("Error: unsupported neighborhood name")
                return
        # translate neighborhood into numpy ndarray
        ca.__neighborhood = np.array(neighborhood, dtype=np.intp)
        # number of cells composing a neighborhood
        ca.size       = np.uint(ca.__neighborhood.shape[0])
        # number of dimmensions
        ca.dimensions = ca.__neighborhood.shape[1]
        # weights for conversion between neighborhood and cell
        ca.weights    = ca.stateset ** np.flip(np.arange(ca.size))
        # neighborhood index array
        ca.index      = np.swapaxes(ca.__neighborhood,0,1)
        # size of the overlap set
        ca.overlapset = ca.stateset**(ca.size-np.uint(1))
        # size of subset diagram
        ca.subset     = np.intp(2)**ca.stateset

    neighborhood = property(get_neighborhood, set_neighborhood)

    def get_rule(ca):
        return ca.__rule

    def set_rule(ca, rule):
        ca.__rule = rule
        # build transition lookup table
        ca.transition_table = np.array([ca.rule // ca.stateset**n % ca.stateset for n in np.arange(ca.stateset**ca.size, dtype=np.uint)], dtype=ca.dtype)
        # for 1D CA some rule properties can be computed
        if (ca.dimensions == 1):
            ca.de_Bruijn = [np.matrix(np.zeros ((ca.overlapset, ca.overlapset), np.uint)) for k in range(ca.stateset)]
            for n in np.arange(ca.stateset**ca.size, dtype=np.uint):
                o_l = n // ca.stateset
                o_r = n % (ca.overlapset)
                ca.de_Bruijn [ca.transition_table[n]] [o_l, o_r] = 1
        # construct the subset diagram
       #ca.Sf = [ [ common.list2int (common.list2bool (np.array(np.mat(common.int2list(i, ca.stateset, ca.overlapset)) * ca.de_Bruijn[c])[0]), ca.stateset) for i in range(ca.subset) ] for c in range(ca.stateset) ]
       #ca.Sb = [ [ common.list2int (common.list2bool (np.array(ca.de_Bruijn[c] * np.mat(common.int2list(i, ca.stateset, ca.overlapset)))[0]), ca.stateset) for i in range(ca.subset)) ] for c in range(ca.stateset) ]

    rule = property(get_rule, set_rule)

    def neighborhood_to_number(ca, neighborhood_array):
        return np.sum(neighborhood_array * ca.weights)

    def neighborhood_to_array (ca, neighborhood_number):
        neighborhood_array = np.zeros(ca.size, dtype=ca.dtype)
        for i in np.flip(np.arange(np.intp(ca.size))):
            neighborhood_array[i] = neighborhood_number % ca.stateset
            neighborhood_number = neighborhood_number / ca.stateset
        return neighborhood_array
    
    def overlap_to_number(ca, overlap_array):
        return np.sum(overlap_array * ca.weights[0:-1])

    def overlap_to_array (ca, overlap_number):
        overlap_array = np.zeros(ca.size-np.uint(1), dtype=ca.dtype)
        for i in np.flip(np.arange(np.intp(ca.size-np.uint(1)))):
            overlap_array[i] = overlap_number % ca.stateset
            overlap_number = overlap_number / ca.stateset
        return overlap_array

    def GoE_count(ca, N):
        M = np.zeros ((ca.subset, ca.subset), int)
        for c in range(ca.stateset):
            for i in range(ca.subset):
                 M[i][ca.Sf[c][i]] = M[i][ca.Sf[c][i]] + 1
        M = np.matrix (M)
        V = np.zeros (ca.subset, int)
        V[ca.subset-1] = 1
        V = np.matrix (V)
        L = []
        for _ in range(N) :
            L.append (V[0,0])
            V = V*M
        #Lf = [float(L[i]/2**i) for i in range(len(L))]  # float GoE density
        return L


class lattice():
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
    1D -- np.ndarray, list or string
    2D -- np.ndarray, list of lists or strings, or from images (png, bmp, ...)
    nD -- np.ndarray only

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
    def __init__(lt, ca, configuration, boundary='cyclic'):
        lt.ca    = ca
        # detecting the CA type
        if (lt.ca.dimensions==1):
            lt.type = '1D'
        else :
            lt.type = 'general'
        # parsing the input configuration
        if (lt.ca.dimensions == 1):
            # 1D CA can be initialized from np.ndarray, list or string
            if (type(configuration) in (np.ndarray, list)):
                lt.configuration = np.array(configuration, dtype=ca.dtype)
            elif (type(configuration) in (str,)):
                lt.configuration = np.fromstring(configuration, dtype=ca.dtype) - ca.dtype(48)
            else:
                print("Error: incompattible 1D configuration")
                return
        elif (lt.ca.dimensions == 2):
            # 2D CA can be initialized from np.ndarray, list of lists or strings, or from images
            if (type(configuration) in (np.ndarray,)):
                lt.configuration = np.array(C, dtype=ca.dtype)
            elif (type(configuration) in (list,)):
                if (type(configuration[0]) in (list,)):
                    lt.configuration = np.array(configuration, dtype=ca.dtype)
                elif (type(configuration) in (string,)):
                    lt.configuration = np.empty(len(configuration), dtype=ca.type)
                    for y in range(len(configuration)):
                        lt.configuration = np.fromstring(configuration[y], dtype=ca.dtype) - ca.dtype(48)
                else:
                    print("Error: incompattible 2D configuration")
                    return
            else:
                print("Error: incompattible 2D configuration")
                return
        else:
            # nD CA where n>2 can only be initialized from np.ndarray
            if (type(configuration) in (np.ndarray,)):
                lt.configuration = np.array(configuration, dtype=ca.dtype)
            else:
                print("Error: incompattible nD configuration")
                return
        lt.shape = lt.configuration.shape
        # an empty array for neighborhoods, can be used for temporal results
        lt.Cn = np.empty(lt.shape, dtype=ca.dtype)
        # parsing the input boundary
        if (boundary in ('cyclic', 'open')):
            lt.boundary = boundary
        else :
            print("Error: currently only 'cyclic' and 'open' boundaries are supported")
            return

    def __repr__(lt):
        return repr((lt.configuration+int(48)).tostring())

    def next(lt, t=1):
        """
        Performs a step forward in time,
        the result is stored back into the configuration.
        """
        if (lt.boundary == 'cyclic'):
            # calculate temporary neighborhoods array
            neighborhoods = np.zeros(lt.configuration.shape, dtype=np.uint)
            for i in np.ndindex(lt.shape):
#                print(np.mod(lt.ca.index + [[i]], lt.shape))
                neighborhood = lt.configuration[np.mod(lt.ca.index + [[i]], lt.shape[0])][0]
#                print(neighborhood)
                neighborhoods[i] = lt.ca.neighborhood_to_number(neighborhood)
#            print(neighborhoods)
            # calculate next configuration from neighborhoods
            # TODO: do not use a for loop for speed, and modify in place
            for i in np.ndindex(lt.shape):
                lt.configuration[i] = lt.ca.transition_table[neighborhoods[i]]
        else:
            print("ERROR: boundary type not supported.")

    def isGoE(lt):
        if ((lt.ca.dimensions == 1) and (lt.boundary == 'open')) :
            s = lt.ca.subset-1
            for x in range(lt.N) : s = lt.ca.Sf[lt.configuration[x]][s];
            return (s == 0)
        else :
            return("Unsupported boundary")

    def repr_preimage_vector(lt, pv):
        return repr(pv)

    def repr_preimage_vector_array(lt, pva):
        return repr(pva)

    def preimage_vector_array_forward(lt, boundary_vector = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            D = [None] * (N+1)
            if (boundary_vector == None):
                D[0] = np.matrix(np.ones(lt.ca.overlapset, dtype=np.uint))
            else:
                D[0] = boundary_vector
            for x in range(N):
                D[x+1] = D[x] * lt.ca.de_Bruijn[lt.configuration[x]]
            return D
        else:
            print("ERROR: only 1D CA supported")

    def preimage_vector_array_backward(lt, boundary_vector = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            D = [None] * (N+1)
            if (boundary_vector == None):
                D[-1] = np.matrix(np.ones(lt.ca.overlapset, dtype=np.uint))
            else:
                D[-1] = boundary_vector
            for x in range(N, 0, -1):
                D[x-1] = (lt.ca.de_Bruijn[lt.configuration[x-1]] * D[x].T).T
            return D
        else:
            print("ERROR: only 1D CA supported")

    def preimage_vector_array(lt, boundary_vector_left = None, boundary_vector_right = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            Df = lt.preimage_vector_array_forward (boundary_vector = boundary_vector_left )
            Db = lt.preimage_vector_array_backward(boundary_vector = boundary_vector_right)
            D = [None] * (N+1)
            for x in range(N+1):
                # Hadamard product
                D[x] = np.matrix(np.multiply(Df[x], Db[x]))
            return D
        else:
            print("ERROR: only 1D CA supported")


    def preimage_matrix_array_forward(lt, boundary_matrix = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            M = [None] * (N+1)
            if (boundary_matrix == None):
                M[0] = np.matrix(np.identity(int(lt.ca.overlapset), dtype=np.uint))
            else:
                M[0] = np.matrix(boundary_matrix)
            for x in range(N):
                M[x+1] = M[x] * lt.ca.de_Bruijn[lt.configuration[x]]
            return M
        else:
            print("ERROR: only 1D CA supported")

    def preimage_matrix_array_backward(lt, boundary_matrix = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            M = [None] * (N+1)
            if (boundary_matrix == None):
                M[-1] = np.matrix(np.identity(int(lt.ca.overlapset), dtype=np.uint))
            else:
                M[-1] = np.matrix(boundary_matrix)
            for x in range(N, 0, -1):
                M[x-1] = lt.ca.de_Bruijn[lt.configuration[x-1]] * M[x]
            return M
        else:
            print("ERROR: only 1D CA supported")

    def preimage_matrix_array(lt, boundary_matrix_left = None, boundary_matrix_right = None):
        if (lt.ca.dimensions == 1):
            N = lt.shape[0]
            Mf = lt.preimage_matrix_array_forward (boundary_matrix = boundary_matrix_left )
            Mb = lt.preimage_matrix_array_backward(boundary_matrix = boundary_matrix_right)
            M = [None] * (N+1)
            for x in range(N+1):
                # Hadamard product
                # TODO: should the result be just a vector?
                M[x] = np.matrix(np.multiply(Mf[x], Mb[x].T))
            return M
        else:
            print("ERROR: only 1D CA supported")


    def previous(lt, boundary_left = None, boundary_right = None):
        if (lt.ca.dimensions == 1):
            N = np.uint(lt.shape[0])

            if (lt.boundary == 'cyclic'):
                # left boundary matrix
                if (boundary_left == None):
                    Ml = np.matrix(np.identity(int(lt.ca.overlapset), dtype=np.uint))
                else:
                    Ml = np.matrix(boundary_left)
                # matrix array
                Mb = lt.preimage_matrix_array_backward(boundary_right)
                # number of preimages
                p = (Ml*Mb[0]).trace()[0,0]
                # array of preimages
                preimages = np.empty((p, lt.shape[0]), dtype=lt.ca.dtype);
                # initialize left side overlaps for each preimage
                o_p0 = np.arange(lt.ca.overlapset, dtype=np.uint).repeat(Mb[0].diagonal().A1.astype(np.intp))
                o_p = o_p0.copy()
               
                for x in np.arange(N, dtype=np.uint):
                    i = 0
                    while (i<p):
                        o_L = o_p[i]; o_R = o_p0[i]
                        for cell in np.arange(lt.ca.stateset, dtype=np.uint):
                            n = o_L * lt.ca.stateset + cell
                            if (lt.configuration[x] == lt.ca.transition_table[n]):
                                o_x = n % lt.ca.overlapset
                                p_i = Mb[x+np.uint(1)][o_x, o_R]
                                for _ in range(p_i):
                                    preimages[i, np.mod(x+np.uint(1), N)] = cell
                                    o_p[i] = o_x
                                    i = i+1

                return preimages

            elif (lt.boundary == 'open'):
                # left boundary vector
                if (boundary_left == None):
                    Dl = np.matrix(np.ones(lt.ca.overlapset, dtype=np.uint))
                else:
                    Dl = np.matrix(boundary_vector)
                # vector array
                Db = lt.preimage_vector_array_backward(boundary_right)
                # number of preimages
                p = (Dl*Db[0].T).sum()
                # array of preimages
                preimages = np.empty((p, lt.shape[0]+int(lt.ca.size)-1), dtype=lt.ca.dtype);
                # left side overlaps for each preimage
                o_p = np.array(np.arange(lt.ca.overlapset), dtype=np.uint).repeat(Db[0].A1.astype(np.intp))
                # initialize left side overlaps for each preimage
                for i in range(p) :
                    preimages[i][0:np.intp(lt.ca.size-np.uint(1))] = lt.ca.overlap_to_array(o_p[i])

                for x in np.arange(N, dtype=np.uint):
                    i = 0
                    while (i<p):
                        o_L = o_p[i];
                        for cell in np.arange(lt.ca.stateset, dtype=np.uint):
                            n = o_L * lt.ca.stateset + cell
                            if (lt.configuration[x] == lt.ca.transition_table[n]) :
                                o_x = n % (lt.ca.overlapset)
                                p_i = Db[x+np.uint(1)].A1[o_x]
                                for _ in range(p_i) :
                                    preimages[i][x+lt.ca.size-+np.uint(1)] = cell
                                    o_p[i] = o_x
                                    i = i+1

                return preimages

            else :
                print("Error: there is no preimage implementation for automata with more than 1 dimensions.")