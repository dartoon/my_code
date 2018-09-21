from pylab import *

# Generate some random data that looks like yours
N = 1000
X = random(N)
Y = sin(X*5) + X*random(N)*.8
Z = random(N)
ERR = X*random(N)

# These are the new arguments that I used
scatter_kwargs = {"zorder":100}
error_kwargs = {"lw":.5, "zorder":0}

scatter(X,Y,c=Z,**scatter_kwargs)
errorbar(X,Y,yerr=ERR,fmt=None, marker=None, mew=0,**error_kwargs )
xlim(0,1)
show()
