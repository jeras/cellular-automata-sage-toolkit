Creation of a CA rule object:

```
sage: import pyca
sage: ca = pyca.rule (2, (-1,0,1), 110)
sage: ca
Cellular Automaton (states = 2, neighborhood = ((-1,), (0,), (1,)), rule = 110)
sage: ca.f
array([0, 1, 1, 1, 0, 1, 1, 0], dtype=uint8)
```

Creation of 1D lattice, three different data types are used for the
configuration, the results are checked to be equal:
```
sage: lt = pyca.lattice (ca, [0,0,0,1,0,0,1,1,0,1,1,1,1,1])
sage: lt
'00010011011111'
sage: for _ in xrange(7) : lt.next(); lt
'00110111110001'
'01111100010011'
'11000100110111'
'01001101111100'
'11011111000100'
'11110001001101'
'00010011011111'
```

Creation of an image of a 1D CA run. The initial configuration is 1024
cells with random value, the preciously defined rule 110 CA is used.
The output image is stored into a local file "rule110_rqndom.png"

Test of neighborhood input format
```
sage: ca_0 = pyca.rule (2, (-1,0,1), 110)
sage: ca_1 = pyca.rule (2, [-1,0,1], 110)
sage: ca_2 = pyca.rule (2, [[-1],[0],[1]], 110)
sage: ca_3 = pyca.rule (2, 'elementary', 110)
sage: ca_0.a == ca_1.a == ca_2.a == ca_3.a
True
```

Changing the neighborhood and rule property of a rule object
```
sage: ca = pyca.rule (2, 'elementary', 110)
sage: ca
Cellular Automaton (states = 2, neighborhood = ((-1,), (0,), (1,)), rule = 110)
sage: ca.a = (1,0,-1)
sage: ca
Cellular Automaton (states = 2, neighborhood = ((1,), (0,), (-1,)), rule = 110)
sage: ca.f
array([0, 1, 1, 1, 0, 1, 1, 0], dtype=uint8)
sage: ca.r = 30
sage: ca.f
array([0, 1, 1, 1, 1, 0, 0, 0], dtype=uint8)
```
