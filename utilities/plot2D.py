import matplotlib.pyplot as plt
import sys

file = open(sys.argv[1])
data = file.read()
lines = data.split('\n')

xs = []
ys = []
zs = []

for i in range(0,len(lines)):
    line = lines[i]
    splt = line.split(',')
    try:
        x = splt[0].strip()
        y = splt[1].strip()
        z = splt[2].strip()
        xs.append(float(x))
        ys.append(float(y))
        zs.append(float(z))
    except:
        continue

plt.tricontourf(xs,ys,zs,100)
plt.show()