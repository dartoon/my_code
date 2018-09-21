import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pickle
import corner

select= input('of which models results to be presented? 0--LCDM;1--wCDM;2:CPL:   ')
if select == 0:
   ndim, nwalkers =2, 50
   sampler,burn=pickle.load(open('int_LCDM_500.txt','rb'))
   plt.plot(burn[0,50:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
   #plt.axhline(4.294, color='blue')
   samples=sampler[:, 50:, :].reshape((-1, ndim))
   print samples.shape
   fig = corner.corner(samples, labels=["$om$", "$h0$"], truths=[0.3,70],
                       quantiles=[0.16, 0.5, 0.84],  #quantiles shows 1D possiblity region
                       show_titles=True, title_kwargs={"fontsize": 12}) 
if select == 1:
   ndim, nwalkers =3, 50
   sampler,burn=pickle.load(open('int_wCDM_500.txt','rb'))
   plt.plot(burn[0,50:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
   #plt.axhline(4.294, color='blue')
   samples=sampler[:, 50:, :].reshape((-1, ndim))
   print samples.shape
   fig = corner.corner(samples, labels=["$om$", "$w$", "$h0$"], truths=[0.3, -1, 70],
                       quantiles=[0.16, 0.5, 0.84],  #quantiles shows 1D possiblity region
                       show_titles=True, title_kwargs={"fontsize": 12}) 
if select == 2:
   ndim, nwalkers =4, 50
   sampler,burn=pickle.load(open('int_CPL_500.txt','rb'))
   plt.plot(burn[0,:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
   #plt.axhline(4.294, color='blue')
   samples=sampler[:, 50:, :].reshape((-1, ndim))
   print samples.shape
   fig = corner.corner(samples, labels=["$om$", "$w0$", "wa", "$h0$"], truths=[0.3, -1, 0, 70],
                       quantiles=[0.16, 0.5, 0.84],  #quantiles shows 1D possiblity region
                       show_titles=True, title_kwargs={"fontsize": 12})
plt.show()
