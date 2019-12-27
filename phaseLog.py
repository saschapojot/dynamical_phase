from sympy import *
import numpy as np
import time
from multiprocessing import Pool

ak, bk, t = symbols('ak,bk,t', cls=Symbol, real=True)
zk = ak + I * bk
z0 = exp(I * zk * t)
z0 = z0.expand(complex=True)
# pprint(z0)
z0c = exp(-I * zk * t)
z0c = z0c.expand(complex=True)
# pprint(z0c)
# calculation of phase
# params
hx, hz, G = symbols('hx,hz,G', cls=Symbol, real=True)
k = symbols('k', cls=Symbol, real=True)
t = symbols('t', cls=Symbol, real=True)

b1, b2, b3, b4 = symbols('b1,b2,b3,b4', cls=Symbol)
rk, sk = symbols('rk,sk', cls=Symbol, real=True)
# zk=symbols('zk',cls=Symbol)
ak, bk = symbols('ak,bk', cls=Symbol, real=True)
R = sqrt(sqrt((-G ** 2 / 4 + hx ** 2 + hz ** 2) ** 2 + (G * hz) ** 2))
solBase = Matrix([[b1, b3], [b2, b4]])
D = Matrix([[z0, 0], [0, z0c]])
psiki = solBase * D * solBase.inv() * Matrix(2, 1, [rk, sk])
Hfk = Matrix([[hz + I * G / 2, hx], [hx, -(hz + I * G / 2)]])

stheta = simplify(G * hz / R ** 2)
ctheta = simplify((-G ** 2 / 4 + hx ** 2 + hz ** 2) / R ** 2)

# theta in [-pi,pi]
cthetaH = sqrt((1 + ctheta) / 2)
sthetaH = sqrt((1 - ctheta) / 2) * sign(stheta)
# phase1=phase1.subs(ak,R*cthetaH)
# phase1=phase1.subs(bk,R*sthetaH)
##########################
# Gkm
Gkm = (Matrix(1, 2, [rk, sk]) * solBase * D * solBase.inv() * Matrix(2, 1, [rk, sk]))[0]

Gkm = Gkm.subs(b1, -hx)
Gkm = Gkm.subs(b2, I * G / 2 + hz + zk)
Gkm = Gkm.subs(b3, -hx)
Gkm = Gkm.subs(b4, I * G / 2 + hz - zk)
Gkm = Gkm.subs(rk, hx / sqrt(hx ** 2 + (hz - sqrt(hx ** 2 + hz ** 2)) ** 2))
Gkm = Gkm.subs(sk, (sqrt(hx ** 2 + hz ** 2) - hz) / sqrt(hx ** 2 + (sqrt(hx ** 2 + hz ** 2) - hz) ** 2))
Gkm = Gkm.subs(ak, R * cthetaH)
Gkm = Gkm.subs(bk, R * sthetaH)

# case1
u = 0.6
r = 0.3
# G=1

Gkm = Gkm.subs(hx, u + r * cos(k))
Gkm = Gkm.subs(hz, r * sin(k))
Gkm = Gkm.subs(G, 1)
Gkm = lambdify([k, t], Gkm, 'numpy')

kStart = 0
kEnd = 2 * np.pi
kN = 400 + 1
kDelta = (kEnd - kStart) / ((kN - 1))
kRange = range(0, kN)
kStep = (kEnd - kStart) / (kN - 1)

tStart = 0
tN = 400 + 1
tStep = 0.05
tRange = range(0, tN)
tDelta = tStep


def f(x):
    xTmp = np.angle(x) / (2 * np.pi)
    if xTmp >= 1 / 2:
        xTmp -= 1
    if xTmp <= -1 / 2:
        xTmp += 1
    return xTmp


def g(x):
    return np.angle(x) / (2 * np.pi)


def phaseLog(argTK):
    tTmp = tStart + argTK[0] * tStep
    kTmp = kStart + argTK[1] * kStep
    valTmp = g(Gkm(kTmp + kStep, tTmp)) - g(Gkm(kTmp, tTmp))
    return [argTK[0], argTK[1], valTmp]


tkArgs = [(j1, j2) for j1 in tRange
          for j2 in kRange]
t1=time.time()
pool=Pool(8)
rL=pool.map(phaseLog,tkArgs)
pool.close()
pool.join()
t2=time.time()
print('phase1 par time: {}'.format(t2-t1))
np.savetxt('jump1.txt',rL)

