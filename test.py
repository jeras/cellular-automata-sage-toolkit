import pyca

ca = pyca.rule (2, (-1,0,1), 110)
print(ca)
lt = pyca.lattice (ca, [0,0,0,1,0,0,1,1,0,1,1,1,1,1])
print(lt)
for _ in range(7) :
    lt.next()
    print(lt)
