import matplotlib.pyplot as plt
import sys

n=1
file = open(sys.argv[n])
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
        xs.append(float(x)*1E10)
        ys.append(float(y)*1E10)
    except:
        continue


plt.hist2d(xs,ys,bins=100)

plt.axes().set_aspect('equal', 'datalim')
plt.xlabel("x (Å)", fontsize=20)
plt.ylabel("z (Å)", fontsize=20)
#plt.xlim(-2,2)
plt.ylim(-3,3)
plt.tight_layout()
plt.show()
