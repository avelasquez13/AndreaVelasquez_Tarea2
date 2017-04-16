import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = np.loadtxt('sedov.dat')
tiempos = np.loadtxt('tiempo.dat')

r, rho_10, rho_60, rho_120 = data[0], data[1], data[2], data[3]
t_10, t_60, t_120 = tiempos[0], tiempos[1], tiempos[2]

i=0
while i<len(r):
    if rho_10[i] == -1:
        r = np.delete(r, i)
        rho_10 = np.delete(rho_10, i)
        rho_60 = np.delete(rho_60, i)
        rho_120 = np.delete(rho_120, i)
	i-=1
    i+=1


plt.plot(r, rho_10, label = r'$t=%d\mathrm{s}$'%t_10, color='yellow')
plt.plot(r, rho_60, label = r'$t=%d\mathrm{s}$'%t_60, color='orange')
plt.plot(r, rho_120, label = r'$t=%d\mathrm{s}$'%t_120, color='r')
plt.xlabel(r'$r\ \mathrm{(m)}$')
plt.ylabel(r'$\rho\ $'+r'$\mathrm{(\frac{kg}{m^3})}$')
plt.title(r'$\mathrm{Densidad\ radial\ promedio}$')
plt.legend()
plt.savefig('prom_rho.png')
plt.close()
