import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astropy.table import Table

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def percentile16(x):
    return np.percentile(x, 16)

def percentile84(x):
    return np.percentile(x, 84)

d = Table.read('HSC_sdss_li21.ascii', format='ascii')


import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'


fig = plt.figure(figsize=(11, 7))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

# HSC: rest 5000 A
size, z, _ = stats.binned_statistic(d['z'], d['rkpc'], bins=np.arange(0.15, 1.0, 0.15),
                                    statistic='median')
ax1.plot(z[:-1], size, 'k', label='HSC (Li+21)')

size_16, z, _ = stats.binned_statistic(d['z'], d['rkpc'], bins=np.arange(0.15, 1.0, 0.15),
                                    statistic=percentile16)
size_84, z, _ = stats.binned_statistic(d['z'], d['rkpc'], bins=np.arange(0.15, 1.0, 0.15),
                                    statistic=percentile84)
ax1.fill_between(z[:-1], size_16, size_84, alpha=0.25)

# JWST ceers
ax1.scatter(1.646, 4.05, color='tab:orange', marker='v') #F150
ax1.scatter(2.317, 1.45, color='tab:orange', marker='v') #F160
ax1.scatter(2.588, 1.11, color='tab:orange', marker='v') #F160
ax1.scatter(3.442, 1.18, color='tab:orange', marker='v') #F200
ax1.scatter(3.465, 1.17, color='tab:orange', label='CEERS (Ding+22)', marker='v') #F200

# JWST SHELLQs: F356W
ax1.scatter(6.34, 1.87, label='This Work', marker='*', s=200, color='r')
ax1.scatter(6.40, 0.75, marker='*', s=200, color='r')

# HST: rest 5000 A
hst_size = np.array([1.17873376, 0.85913774, 3.62829472, 1.95819451, 4.36564811,
       3.32197847, 2.70039055, 2.18717427, 7.55972388, 1.76052273,
       6.05772044, 3.12287316, 3.42177231, 1.4186286 , 1.77963679,
       2.62772947, 1.55200827, 5.44468214, 0.87061829, 1.63473507,
       2.27353186, 5.87892781, 0.95677881, 4.33088298, 0.881597  ,
       5.066965  , 2.50065858, 3.10583291, 0.99854718, 1.48723315,
       2.06544458, 1.63706414])

hst_z = np.array([1.63 , 1.301, 1.667, 1.447, 1.326, 1.57 , 1.552, 1.567, 1.618,
       1.532, 1.244, 1.407, 1.478, 1.239, 1.294, 1.617, 1.527, 1.579,
       1.325, 1.411, 1.276, 1.412, 1.585, 1.551, 1.516, 1.6  , 1.483,
       1.626, 1.337, 1.272, 1.445, 1.664])

ax1.scatter(hst_z, hst_size, color='tab:cyan', label='HST (Ding+20)')
ax1.set_xlabel('$z$', fontsize=30)
ax1.set_ylabel(u'$R_{\mathrm{eff, maj}}$ (kpc)', fontsize=30)
# u'log($R_{\mathrm{eff, circ.}}$, units of kpc)

ax1.tick_params(labelsize=20)
# BlueTides: Sersic radius, F356W
blue = Table.read('bluetides.txt', format='ascii')
ax1.scatter([7.0]*len(blue), blue['size_sersic'], label='BlueTides (Marshall+19)', marker='^', color='gray', facecolor='None')
ax1.set_xlim(-0.1, 7.5)

# galaxy
# Mowla+19: rest 5000 A
mowla_z = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
sf = 6.3*(1+mowla_z)**(-0.5)  # 10.25 < logM* < 10.75
qs = 4.1*(1+mowla_z)**(-1.1)

ax1.plot(mowla_z, sf, label='SF galaxy (Mowla+19)', color='b', linestyle='dashed')
ax1.plot(mowla_z, qs, label='QS galaxy (Mowla+19)', color='r', linestyle='dashdot')

plt.tick_params(labelsize=20)

# Allen+17: F160W
allen_z = np.arange(1, 7.1, 0.1)
sf = 7.07*(1+allen_z)**(-0.89)  # logM* > 10
ax1.plot(allen_z, sf, label='SF galaxy (Allen+17)', color='gray', linestyle='dotted')

# look back time
lookback = cosmo.lookback_time(np.arange(0, 8, 2))
lookback = ['{:.1f}'.format(time.value) for time in lookback]

new_tick_locations = np.arange(0, 8, 2)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(lookback)
ax2.set_xlabel('lookback time (Gyr)', fontsize=25)

plt.tick_params(labelsize=20)

ax1.legend(frameon=True, fontsize=16)

fig.savefig('quasar_host_size.pdf')
