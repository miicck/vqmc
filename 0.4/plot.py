import matplotlib.pyplot as plt
import sys
file = open(sys.argv[1])
allData = file.read()
datas = allData.split('#')

xss = []
yss = []

for d in range(0,len(datas)):

    data = datas[d]
    split = data.split('\n')

    xs = []
    ys = []

    for i in range(0,len(split)):
        entry = split[i]
        pair = entry.split(',')
        try:
            x = pair[0].strip()
            y = pair[1].strip()
            xs.append(float(x))
            ys.append(float(y))
        except:
            continue
    
    xss.append(xs)
    yss.append(ys)

for d in range(0,len(xss)):
    plt.plot(xss[d],yss[d])
#plt.ylim(-3,10)
#plt.xlabel("c1 = 1-c2")
#plt.ylabel("Energy (eV)")
plt.axes().set_aspect('equal', 'datalim')
plt.show()