from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["font.family"] = "sans-serif"
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(111)

d = Table.read('prob.txt', format='ascii')

prob = d['prob']
gamma = d['gamma']

dx = gamma[1]-gamma[0]
prob_sum = np.sum(prob * dx)

ax.plot(gamma, prob)
ax.set_xlabel('$\gamma$', fontsize=12)
ax.set_ylabel('pdf', fontsize=12)

ax.axvline(0.22, 0, 1, color='k', linestyle='dashed')
ax.fill_between([0.22-0.62, 0.22+0.63], 0, 1, color='silver', alpha=0.3, edgecolor='None')
ax.set_ylim(0, 1)

ax.tick_params(axis='x', which='both', direction='in', left=True, bottom=True, right=True, top=True, length=2.5)
ax.tick_params(axis='y', which='both', direction='in', left=True, bottom=True, right=True, top=True, length=2.5)

fig.subplots_adjust(top=0.95, right=0.98, left=0.13, bottom=0.14)
fig.savefig('gamma.pdf')
