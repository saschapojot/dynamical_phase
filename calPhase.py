import numpy as np
import matplotlib.pyplot as plt

kStart = 0
kEnd = 2 * np.pi
kN = 400 + 1
kRange = range(0, kN)
kStep = (kEnd - kStart) / (kN - 1)

tStart = 0
tN = 400 + 1
tStep = 0.05
tRange = range(0, tN)
jump1 = np.array(np.loadtxt('jump1.txt'))
intPhase = np.array(np.loadtxt('phaseCorrected.txt'))

jump1 = jump1[jump1[:, 1].argsort(kind='mergesort')]
jump1 = jump1[jump1[:, 0].argsort(kind='mergesort')]

intPhase = intPhase[intPhase[:, 1].argsort(kind='mergesort')]
intPhase = intPhase[intPhase[:, 0].argsort(kind='mergesort')]








intPhaseSortByK=intPhase[intPhase[:,1].argsort(kind='mergesort')]
b=[]
C=intPhaseSortByK[:,2]
seulK=np.unique(intPhaseSortByK[:,1])
for kElem in seulK:
    thisCols=C[np.where(intPhaseSortByK[:,1]==kElem)]
    for j in range(0,tN):
        if j==0:
            b.append(0)
        else:
            b.append(sum(thisCols[0:j]))

preB=np.zeros((len(intPhaseSortByK),len(intPhaseSortByK[0])))
preB[:,0]=intPhaseSortByK[:,0]
preB[:,1]=intPhaseSortByK[:,1]
preB[:,2]=b
preB=preB[preB[:,0].argsort(kind='mergesort')]
nRowPreB=len(preB)
deltaB=[]
zeit=np.unique(intPhase[:,0])
for tElem in zeit:
    thisCols=preB[np.where(preB[:,0]==tElem)][:,2]
    thisLen=len(thisCols)

    thisDelta=thisCols[1:(thisLen)]-thisCols[0:(thisLen-1)]
    deltaB=np.concatenate([deltaB,thisDelta])



B=np.zeros((((kN-1)*tN),3))
tkRows=[(j1,j2) for j1 in range(0,tN) for j2 in range(0,kN-1)]
for j3 in range(0,len(tkRows)):
    B[j3,0]=tkRows[j3][0]
    B[j3,1]=tkRows[j3][1]
    B[j3,2]=deltaB[j3]


deltaA=[]
for tElem in zeit:
    thisCols=jump1[np.where(jump1[:,0]==tElem)][:,2]
    thisLen=len(thisCols)
    thisDelta=thisCols[1:thisLen]-thisCols[0:(thisLen-1)]
    deltaA=np.concatenate([deltaA,thisDelta])
A=np.zeros((((kN-1)*tN),3))
for j3 in range(0,len(tkRows)):
    A[j3,0]=tkRows[j3][0]
    A[j3,1]=tkRows[j3][1]
    A[j3,2]=deltaA[j3]
threshold=0.1
def f(x):
    if x >= threshold:
        x -= 1
    if x <= -threshold:
        x += 1
    return x
vdelta=np.zeros((((kN-1)*tN),3))
for j3 in range(0,len(tkRows)):
    vdelta[j3,0]=tkRows[j3][0]
    vdelta[j3,1]=tkRows[j3][1]
    vdelta[j3,2]=f(deltaA[j3]+deltaB[j3])

vGeo=[]
for tElem in zeit:
    vGeo.append(sum(vdelta[np.where(vdelta[:,0]==tElem)][:,2]))


plt.plot(np.array((tRange))*0.05,vGeo)
plt.show()



