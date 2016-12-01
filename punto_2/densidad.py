import numpy as np
import corner
import matplotlib.pyplot as plt

log10M=np.genfromtxt("log10M.txt")
alpha=np.genfromtxt("alpha.txt")

ll=[]
aa=[]

for i in range(len(log10M)):
    li=[log10M[i]]
    ai=[alpha[i]]
    ll.append(li)
    aa.append(ai)

data=np.hstack((ll,aa))

figure = corner.corner(data, labels=[r"$log_{10}M_{sol}$", r"$\alpha$", r"$\Gamma \, [\mathrm{parsec}]$"],quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
filep='densidad.pdf'
plt.savefig(filep,format='pdf')
plt.close()
