import matplotlib.pyplot as plt


file = open("log")
data = file.read()
splt = data.split(',')
pts = []

for i in range(0,len(splt)):
    st = splt[i]
    try:
        val = float(st)
        pts.append(val)
    except:
        continue

plt.hist(pts,bins=500)
plt.yscale('log')
plt.show()