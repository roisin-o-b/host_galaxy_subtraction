import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import sys

counts = ascii.read('counts.csv') #read in STELLA data
obs = Table(names = ('dec_year', 'av1', 'err1', 'av2', 'err2', 'av3', 'err3', 'avmrk1018', 'errmrk1018')) #create new table for averaged exposures

#calculate error-weighted averages
dates = np.unique(counts['dec_year'])

for d in dates:
    bln = counts['dec_year'] == d #boolean for all exposures on the same day
    av1 = np.average(counts[bln]['ref1'], weights= 1/counts[bln]['err1']**2) # error weighted exp averages for stars and mrk1018
    av2 = np.average(counts[bln]['ref2'], weights= 1/counts[bln]['err2']**2)
    av3 = np.average(counts[bln]['ref3'], weights= 1/counts[bln]['err3']**2)
    avA = np.average(counts[bln]['mrk1018' ], weights= 1/counts[bln]['err_mrk1018']**2)
    
    err1 = np.sqrt(1/np.sum(1/counts[bln]['err1']**2)) #errors for weighted averages
    err2 = np.sqrt(1/np.sum(1/counts[bln]['err2']**2))
    err3 = np.sqrt(1/np.sum(1/counts[bln]['err3']**2))
    err_mrk1018 = np.sqrt(1/np.sum(1/counts[bln]['err_mrk1018']**2))
    
    obs.add_row([d, av1, err1, av2, err2, av3, err3, avA, err_mrk1018])
    
ascii.write(obs, 'new_obs.csv', format = 'csv', overwrite = True)

#------------------------plot differential LCs-------------------------------

zp = 15.87

#convert counts to mags
mags = Table(data = (-2.5*np.log10(obs['av1']),  -2.5*np.log10(obs['av2']), -2.5*np.log10(obs['av3']), -2.5*np.log10(obs['avmrk1018'])), names = ('r1', 'r2', 'r3', 'mrk1018'))

mags_err = Table(data = (1.09*(obs['err1']/obs['av1']), 1.09*(obs['err2']/obs['av2']), 1.09*(obs['err3']/obs['av3']),1.09*(obs['errmrk1018']/obs['avmrk1018'])), names = ('r1', 'r2', 'r3', 'mrk1018'))

mags_corr = Table(data = (obs['dec_year'], mags['r1'] - mags['r3']+zp, mags['r2'] - mags['r3']+zp, mags['r3'] - mags['r3']+zp, mags['mrk1018'] - mags['r3']+zp), names = ('dec_year', 'r1_corr', 'r2_corr', 'r3_corr', 'mrk1018_corr'))

mags_corerr = Table(data=(np.sqrt(mags_err['r1']**2+mags_err['r3']**2), np.sqrt(mags_err['r2']**2+mags_err['r3']**2), np.sqrt(mags_err['r3']**2+mags_err['r3']**2), np.sqrt(mags_err['mrk1018']**2+mags_err['r3']**2)), names = ('r1', 'r2', 'r3', 'mrk1018'))

ascii.write(mags_corerr, 'mrk1018_magser.csv', format = 'csv', overwrite = True)
ascii.write(mags_corr, 'mrk1018_mags.csv', format = 'csv', overwrite = True)

#--------------------------correct for zeropoint mag----------------------


plt.gca().invert_yaxis()
#plt.ylim(17.5, 16.3)
plt.ylabel('u-band magnitude')
plt.xlabel('decimal year')
plt.errorbar(obs['dec_year'], mags_corr['mrk1018_corr'], yerr = mags_corerr['mrk1018'], linestyle='None', marker='.', color='purple')
plt.title('Mrk 1018 AGN + host lightcurve')
plt.ticklabel_format(useOffset=False)
plt.show()
plt.savefig('combined_lightcurve')
plt.close()
