import numpy as np
import matplotlib.pyplot as plt



xh=np.arange(0,24,1)
xg=np.roll(xh,-15)
mk=np.roll(xh,-9)
xg=xg[(xh<1)+(xh>7)]
mk=mk[(xh<1)+(xh>7)]
xh=xh[(xh<1)+(xh>7)]
print 'xh time in L.A. :  {0}\r;'.format(xh)
print 'mk time in Porland:{0}\r;'.format(mk)
print 'xg time in BJ:     {0}\r;'.format(xg)
l=len(xh)
plt.plot(xh,xh, 'o')
plt.plot(mk,xh, '*')
plt.plot(xg,xh, 'r.')
plt.plot(np.linspace(0,24,50),8+np.zeros(50))
plt.xlim([0,24])
plt.show()
