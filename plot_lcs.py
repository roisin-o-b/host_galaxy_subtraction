import glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

g_1s=ascii.read('data_points_all_gmos_1sig.csv')
v_1s=ascii.read('data_points_all_vimos_1sig.csv')
pts =ascii.read('data_points_all_vimos_mean.csv')

fig,ax = plt.subplots(1)
ax.set_ylabel('Host-subtracted AGN $u\'$-band magnitude')
ax.set_xlabel('Time, decimal year')
ax.set_ylim(22.75, 17)

for i in range(len(pts)):
    if np.isnan(pts[i]['agn_abs']):
        pass
    elif pts[i]['lower_err'] != 0:
        ax.errorbar(pts[i]['dec_yr'], pts[i]['agn_abs'], yerr=[[pts[i]['lower_err']], [pts[i]['upper_err']]], fmt='.b', label = 'AGN only')
    else: 
        lolims = np.ones(1, dtype=bool)
        ax.errorbar(pts[i]['dec_yr'], pts[i]['agn_abs'], yerr=0.3, lolims = lolims, fmt='_b')

for i in range(len(g_1s)):
    if g_1s[i]['agn_abs'] > pts[i]['agn_abs']:
        g_1s[i]['agn_abs'] = pts[i]['agn_abs']
    if np.isnan(g_1s[i]['agn_abs']):
        g_1s[i]['agn_abs'] = pts[i]['agn_abs']
    if np.isnan(v_1s[i]['agn_abs']):
        v_1s[i]['agn_abs'] = 100
    if g_1s[i]['lower_err'] == 0:
        v_1s[i]['agn_abs'] = 100
    if g_1s[i]['lower_err'] == 0 and pts[i]['agn_abs'] > g_1s[i]['lower_err']:
        g_1s[i]['agn_abs'] = pts[i]['agn_abs']

ax.fill_between(g_1s['dec_yr'][0:7], g_1s['agn_abs'][0:7], v_1s['agn_abs'][0:7], color = 'blue', alpha = 0.15)
ax.fill_between(g_1s['dec_yr'][7:], g_1s['agn_abs'][7:], v_1s['agn_abs'][7:], color = 'blue', alpha = 0.15)

#sunblock period
ax.axvspan(2020.2, 2020.4, color='grey', alpha = 0.5)
ax.text(2020.27, 21, 'sunblock', size = 17, rotation = 90)
ax.axvspan(2021.2, 2021.4, color='grey', alpha=0.5)
ax.text(2021.27, 21, 'sunblock', size = 17, rotation = 90)

plt.savefig('lc_paper_2sunbl')
plt.show()
