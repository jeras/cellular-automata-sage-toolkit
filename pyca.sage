def int2list (number, radix, length) :
  list = Integer(number).digits(radix)
  return list + (length-len(list))*[0]

def list2int (list, radix) :
  return sum ( [radix^i*list[i] for i in xrange(len(list))] )

def list2bool (list) :
  return [Integer(list[i]>0) for i in xrange(len(list))]


class CA1D_rule () :
  def __init__(ca, k, m, r):
    ca.k = k; ca.m = m; ca.r = r; 
    ca.f = int2list (ca.r, ca.k, ca.k^ca.m)
    ca.D = [zero_matrix(ZZ, ca.k^(ca.m-1), sparse=True) for k in xrange(ca.k)]
    for n in xrange(ca.k^ca.m) :
      o_l = n // ca.k; o_r = n % (ca.k^(ca.m-1))
      ca.D [ca.f[n]] [o_l, o_r] = 1
    ca.Sf = [ [ list2int (list2bool (vector(int2list(i, ca.k, ca.k^(ca.m-1))) * ca.D[c]), \
                ca.k) for i in xrange(2^(ca.k^(ca.m-1))) ] for c in xrange(ca.k) ]

  def __repr__(ca):
    return "Cellular Automaton:\n"+     \
      "  states    = "+str(ca.k)+"\n" + \
      "  neighbors = "+str(ca.m)+"\n" + \
      "  rule      = "+str(ca.r)


class CA1D_lattice () :
  def __init__(lt, ca, C, sh=0, b='cyclic') :
    lt.ca = ca; lt.C = C; lt.sh = sh; lt.N = len(lt.C); lt.b = b

  def __repr__(lt):
    return repr(lt.C)

  def next (lt) :
    if (lt.b == 'cyclic') :
      lt.C = [ lt.ca.f[list2int([lt.C[(x+lt.ca.m-1-i-lt.sh) % lt.N] \
               for i in xrange(lt.ca.m)], lt.ca.k)] for x in xrange(lt.N            ) ]
    else :
      lt.C = [ lt.ca.f[list2int([lt.C[ x+lt.ca.m-1-i              ] \
               for i in xrange(lt.ca.m)], lt.ca.k)] for x in xrange(lt.N-(lt.ca.m-1)) ]
      lt.N = lt.N - (lt.ca.m-1)
    return

  def prev (lt) :
    C_p = []

    if (lt.b == 'cyclic') :
      lt.D_x_b = [identity_matrix(ZZ, lt.ca.k^(lt.ca.m-1), sparse=True)]
      for x in xrange(lt.N) :
        lt.D_x_b.append (lt.ca.D[lt.C[lt.N-1-x]] * lt.D_x_b[x])
      lt.p = sum ([lt.D_x_b [lt.N] [i,i] for i in xrange(lt.ca.k^(lt.ca.m-1))])

      C_p = [CA1D_lattice(lt.ca, lt.N*[0], lt.sh, lt.b) for i in xrange(lt.p)]
      o_p0 = [];
      for o in xrange(lt.ca.k^(lt.ca.m-1)) :
        o_p0.extend([o for d in xrange(lt.D_x_b[lt.N][o,o])])
      o_p = list(o_p0)

      for x in xrange(lt.N) :
        i = 0
        while (i<lt.p) :
          o_L = o_p[i]; o_R = o_p0[i]
          for c in xrange(lt.ca.k) :
            n = o_L * lt.ca.k + c
            if (lt.C[x] == lt.ca.f[n]) :
              o_x = n % (lt.ca.k^(lt.ca.m-1))
              p_i = lt.D_x_b[lt.N-x-1][o_x,o_R]
              for p_c in xrange(p_i) :
                C_p[i].C [(x+lt.sh) % lt.N] = c
                o_p[i] = o_x
                i = i+1

    else :
      if (lt.b == 'open') :
        b_L = b_R = vector(lt.ca.k^(lt.ca.m-1)*[1])
        lt.b_x_b = [b_R]
      else :
        b_L = vector(lt.b[0]); b_R = vector(lt.b[1])
        lt.b_x_b = [b_R]

      for x in xrange(lt.N) :
        lt.b_x_b.append (lt.ca.D[lt.C[lt.N-1-x]] * lt.b_x_b[x])
      lt.p = b_L * lt.b_x_b[lt.N-1]

      C_p = [CA1D_lattice(lt.ca, (lt.N+ca.m-1)*[0], lt.sh, lt.b) for i in xrange(lt.p)]
      o_p = [];
      for o in xrange(lt.ca.k^(lt.ca.m-1)) :
        o_p.extend([o for d in xrange(b_L[o] * lt.b_x_b[lt.N][o])])
      for i in xrange(lt.p) :
        C_p[i].C [0:lt.ca.m-1] = int2list(o_p[i], lt.ca.k, lt.ca.m-1)

      for x in xrange(lt.N) :
        i = 0
        while (i<lt.p) :
          o_L = o_p[i];
          for c in xrange(lt.ca.k) :
            n = o_L * lt.ca.k + c
            if (lt.C[x] == lt.ca.f[n]) :
              o_x = n % (lt.ca.k^(lt.ca.m-1))
              p_i = lt.b_x_b[lt.N-x-1][o_x]
              for p_c in xrange(p_i) :
                C_p[i].C [x+lt.ca.m-1] = c
                o_p[i] = o_x
                i = i+1

    return C_p

  def isGoE (lt) :
    if (lt.b == 'open') :
      s = 2^(ca.k^(ca.m-1))-1
      for x in xrange(lt.N) : s = lt.ca.Sf[lt.C[x]][s]
      return (s == 0)
    else :
      return "Unsupported boundary"
