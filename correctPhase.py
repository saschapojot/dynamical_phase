import numpy as np
par=np.loadtxt('par.txt')
correction=np.loadtxt('correction.txt')
correction[:,0]=correction[:,0].astype(int)
correction[:,1]=correction[:,1].astype(int)
nRowCor=len(correction)
print(nRowCor)
for j in range(0,nRowCor):
    b1=set(np.where(par[:,0]==correction[j,0])[0])
    b2=set(np.where(par[:,1]==correction[j,1])[0])

    rTmp=list(b1.intersection(b2))[0]
    print(rTmp)
    par[rTmp,2]=correction[j,2]
np.savetxt('phaseCorrected.txt',par)

