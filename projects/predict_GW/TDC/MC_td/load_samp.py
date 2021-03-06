import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import pickle
#sampler=pickle.load(open('samp.txt','rb'))
ndim, nwalkers =2, 100
#samples = sampler.chain[:, 50:, :].reshape((-1, ndim))


#plt.plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#plt.axhline(0.3, color='blue')
samples,burn=pickle.load(open('samp.txt','rb'))
#print np.shape(samples[50*ndim*nwalkers:])

#m_true = -0.9594
#b_true = 4.294
#f_true = 0.534


plt.plot(burn[1], '-', color='k', alpha=0.3)
#plt.axhline(4.294, color='blue')

import corner
fig = corner.corner(samples, labels=["$b0$", "$b1$"],
                      truths=[0.3,70])

plt.show()
