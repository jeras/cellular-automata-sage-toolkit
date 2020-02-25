# a couple of auxiliary functions
def int2list (x, r, n) :
    """
    Translates number x into a list l of n digits base r (LSB first list).
    INPUTS:
        x -- the number to translate
        r -- base
        n -- digits
    OUTPUT:
        l -- list of digits
    """
    l = []
    for i in range(n) : l.append(int(x)%int(r)); x=int(x)/int(r)
    return l

def list2int (l, r) :
    """
    Translates a list l of n digits base r into an integer x (LSB first list).
    INPUTS:
        l -- list of digits
        r -- base
    OUTPUT:
        x -- integer
    """
    x = 0
    for i in range(len(l)) : x = x + l[i] * r**i
    return x

def list2bool (l) :
    for i in range(len(l)) : 
        if (l[i] > 0) : l[i] = 1 
    return [int(l[i]>0) for i in range(len(l))]
