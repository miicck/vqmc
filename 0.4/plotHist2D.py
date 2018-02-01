import matplotlib.pyplot as plt
import sys

file = open(sys.argv[1])
data = file.read()
lines = data.split('\n')

xs = []
ys = []

print("Plotting 2D histogram from file: " + sys.argv[1])

for i in range(0,len(lines)):
    line = lines[i]
    splt = line.split(',')
    try:
        x = splt[0].strip()
        y = splt[1].strip()
        xs.append(float(x))
        ys.append(float(y))
    except:
        continue

plt.hist2d(xs,ys,bins=100)
plt.axes().set_aspect('equal', 'datalim')
plt.show()