import matplotlib.pyplot as plt
import sys

f, axes = plt.subplots(ncols=len(sys.argv)-1)

for n in range(1,len(sys.argv)):
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
            xs.append(float(x))
            ys.append(float(y))
        except:
            continue

    axes[n-1].set_aspect('equal', 'datalim')
    axes[n-1].set_title(str(n) + " total points: " + str(len(xs)))
    axes[n-1].hist2d(xs,ys,bins=100)

plt.show()
