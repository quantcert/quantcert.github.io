import numpy as np

#definition of the Pauli matrices and T the tensor reprensenting each four Paulis' I=T[0], X=T[1], Z=T[2], Y=T[3]
X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]]);Y=np.array([[0,-1j],[1j,0]])

T=[np.identity(2),X,Z,Y]

# convert is a standard function to change number n in decimal basis to the same number in basis b
def convert(n,b,verbose=False):
    t=[]
    if n==0:
        return [0]
    while n>0:
        quotient, remainder = n//b,n%b
        if verbose:
            print("%s=%s*%s+%s"%(n,quotient,b,remainder))
        t.append(remainder)
        n = quotient
    #one adds 0 to make sure the table is of size 2 to correspond to a 2-qubit operator
    if np.size(t)==1:
        t.append(0)
    #one reverses the order of the remainders    
    t.reverse()
    return t

#one computes the 15 points corresponding to the 15 nontrivial two qubit operators. An operator T[i]\otimes T[j] is encoded in the point n=4*i+j i.e by a decimal number n which is equal to ij in base 4
points=[None for i in range (16)]
for i in range(4):
    for j in range(4):
        points[i*4+j]=np.kron(T[i],T[j])
points.remove(points[0]);

#one computes all contexts. On considers all triples of points and only keep the triple whose product is \pm I_4. We also keep track of the sign of each context

All_lines=[]
from itertools import combinations
comb=list(combinations([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],3))

for i in comb:
    M=points[i[0]]@points[i[1]]@points[i[2]]
    #print(M)
    if (M==np.identity(4)).all():
        All_lines.append([i,1])
    if (M==-np.identity(4)).all():
        All_lines.append([i,-1])

#when printing one obtains 15 lines, three of them being negative
print(All_lines)
