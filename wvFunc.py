import matplotlib.pyplot as plt
import numpy as np
y=np.loadtxt('intRst.csv')
NStep = 400
tStart=0
tStep = 0.05
stepRange = range(0, NStep + 1)
timeRange = []
for j in stepRange:
    timeC = tStart + j * tStep
    timeRange.append(timeC)
plt.plot(timeRange,y)
plt.show()
