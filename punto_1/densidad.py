import numpy as np
import corner
import matplotlib.pyplot as plt

x=np.genfromtxt("x.txt")
y=np.genfromtxt("y.txt")

xx=[]
yy=[]

for i in range(len(x)):
    xi=[x[i]]
    yi=[y[i]]
    xx.append(xi)
    yy.append(yi)


data=np.hstack((xx,yy))

figure = corner.corner(data, labels=[r"$x$", r"$y$", r"$\Gamma \, [\mathrm{parsec}]$"],quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
filep='densidad.pdf'
plt.savefig(filep,format='pdf')
plt.close()
