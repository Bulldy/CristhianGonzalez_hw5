import numpy as np
import corner
import matplotlib.pyplot as plt


a=np.genfromtxt("a.txt")
b=np.genfromtxt("b.txt")
c=np.genfromtxt("c.txt")
d=np.genfromtxt("d.txt")

aa=[]
bb=[]
cc=[]
dd=[]

for i in range(len(a)):
    ai=[a[i]]
    bi=[b[i]]
    ci=[c[i]]
    di=[d[i]]
    aa.append(ai)
    bb.append(bi)
    cc.append(ci)
    dd.append(di)

data=np.hstack((aa,bb,cc,dd))

figure = corner.corner(data, labels=[r"$\alpha$", r"$\beta$", r"$\gamma$",r"$\delta$", r"$\Gamma \, [\mathrm{parsec}]$"], quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
filep='densidad.pdf'
plt.savefig(filep,format='pdf')
plt.close()


