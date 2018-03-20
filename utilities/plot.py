from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

def model(x,a,b,c,d,e,f,g):
    return f/(x**3) + e/(x**2) + a/x + b + c*x + d*(x**2) + g*(x**3)


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
    plt.plot(xss[d],yss[d],color="black")
    best_vals, covar = curve_fit(model,xss[d],yss[d])
    print (best_vals)
    ymodel = []
    for i in range(0,len(xss[d])):
        ymodel.append(model(xss[d][i],*best_vals))
    plt.plot(xss[d],ymodel,color="red")

#plt.ylim(-2.875,-2.8)
#plt.axes().set_yscale("log")
#plt.xlim(0,10)
plt.xlabel("1s mixing coefficient (c)", fontsize=20)
plt.ylabel("$E_V$ (Hartree)", fontsize=20)
plt.tight_layout()
#plt.axes().set_aspect('equal', 'datalim')
plt.show()