import numpy as np
X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]]);Y=np.array([[0,-1j],[1j,0]])
T=[np.identity(2),X,Z,Y]

### To express n in base b with 3 bits
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
    #on renverse le tableau des restes successifs
    if np.size(t)==2:
        t.append(0)
    if np.size(t)==1:
        t.append(0)
        t.append(0)
    t.reverse()
    return t

### Create the 63 operators of W(5,2). The point numbered n corresponds to the operator T[i]\otimes T[j]\otimes T[k] where ijk is the base 4 decomposition of n
points=[None for i in range (64)]
for i in range(4):
    for j in range(4):
        for k in range(4):
            points[i*4*4+j*4+k]=np.kron(T[i],np.kron(T[j],T[k]))
points.remove(points[0]);

### Create the 315 signed lines of W(5,2)

All_lines=[]
from itertools import combinations
comb=list(combinations(list(range(0,63)),3))
for i in comb:
    M=points[i[0]]@points[i[1]]@points[i[2]]
    if (M==np.identity(8)).all():
        All_lines.append([i,1])
    if (M==-np.identity(8)).all():
        All_lines.append([i,-1])

print(All_lines)

### Quadrics Q_0(5,2)

### Selection of the points of W(5,2) which satisfies Q_0(x)=0
pointsQp=[]
for i in range(4):
    for j in range(4):
        for k in range(4):
            if i==0:
                x1=[0,0]
            if i==1:
                x1=[0,1]
            if i==2:
                x1 = [1, 1]
            if i==3:
                x1=[1,0]
            if j==0:
                x2=[0,0]
            if j==1:
                x2=[0,1]
            if j==2:
                x2 = [1, 1]
            if j==3:
                x2=[1,0]
            if k==0:
                x3=[0,0]
            if k==1:
                x3=[0,1]
            if k==2:
                x3 = [1, 1]
            if k==3:
                x3=[1,0]
            if ((x1[0]*x1[1]+x2[0]*x2[1]+x3[0]*x3[1]) %2==0):
                pointsQp.append([i,j,k])

pointsQp.remove(pointsQp[0])
print(pointsQp)

### Points of the Quadrics as operators

pointsQ=[None for i in range (35)]
for i in range(35):
    pointsQ[i]=np.kron(T[pointsQp[i][0]],np.kron(T[pointsQp[i][1]],T[pointsQp[i][2]]))

### Lines of the quadrics
All_linesQ=[]
combq=list(combinations(list(range(0,35)),3))
for i in combq:
    M=pointsQ[i[0]]@pointsQ[i[1]]@pointsQ[i[2]]
    #Converstion of the permutation in the original labeling of the lines
    m=(pointsQp[i[0]][0]*4**2+pointsQp[i[0]][1]*4+pointsQp[i[0]][2]-1,pointsQp[i[1]][0]*4**2+pointsQp[i[1]][1]*4+pointsQp[i[1]][2]-1,pointsQp[i[2]][0]*4**2+pointsQp[i[2]][1]*4+pointsQp[i[2]][2]-1)
    if (M==np.identity(8)).all():
        All_linesQ.append([m,1])
    if (M==-np.identity(8)).all():
        All_linesQ.append([m,-1])
        
print(All_linesQ)

