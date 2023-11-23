import glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack
import matplotlib.pyplot as plt




def group_dps(files): 
    
    files = sorted(files)
    dates = np.array(ascii.read(files[0])['dec_yr']) #array of dates in decimal years
    np.savetxt('dec_yrs.csv', dates) #save decimal year dates
    
    agn_only = Table(names = dates) #create empty table for AGN count values for each observation
    
    for phile in files:
        data = ascii.read(phile)
        agn_only.add_row(data['av_agnonly']) #isolate AGN counts from each run
        
    for d in dates:
        dp = np.array(agn_only[str(d)]) #link AGN count with decimal year
        dp = np.sort(dp) #sort data points from lowest to highest
        np.savetxt('data_points/data_point'+str(d)+'.csv', dp) #save counts from same observations in separate files




def agnonlystats(dpts, dec_yrs,):
    
    dpts = sorted(dpts) #sort data points by observation
    agnonlystats = Table(names = ('dec_yr','lower_3sig', 'lower_2sig', 'lower_1sig','centre', 'upper_1sig', 'upper_2sig', 'upper_3sig')) #create table for agn count statsistics
    
    for i, dpt in enumerate(dpts):
        data = np.loadtxt(dpt) #load count values
        dec_yr = dec_yrs[i] #find the obs date
        centre = (data[500]+data[499])/2 #take the central value
        
        upper_3sig = data[998] #upper and lower 3 sigma
        lower_3sig = data[1]
        
        upper_2sig = data[977] #upper and lower 2 sigma
        lower_2sig = data[22]
        
        upper_1sig = data[841] #upper and lower 1 sigma
        lower_1sig = data[158]
        
        agnonlystats.add_row([dec_yr, lower_3sig, lower_2sig, lower_1sig, centre, upper_1sig, upper_2sig, upper_3sig]) #add to table
        
    ascii.write(agnonlystats, 'agnonlystats.csv', format='csv', overwrite = True) #save results




def plot_agnstats(obs, stats):
    
    #find magnitudes of combined host+AGN
    
    zp = 15.87 #reference 3 star absolute magnitude from SDSS
    
    #convert counts to mags from STELLA obs for ref 3 and Mrk 1018 with errors
    
    mags = Table(data = (-2.5*np.log10(obs['av3']), -2.5*np.log10(obs['avmrk1018'])), names = ('r3', 'mrk1018'))
    
    mags_err = Table(data = (1.09*(obs['err3']/obs['av3']),1.09*(obs['errmrk1018']/obs['avmrk1018'])), names = ('r3', 'mrk1018'))
    
    mags_corr = Table(data = (obs['dec_year'], mags['r3'] - mags['r3'], mags['mrk1018'] - mags['r3']), names = ('dec_year', 'r3_corr', 'mrk1018_corr')) #correct for zeropoint
    
    mags_corerr = Table(data=(np.sqrt(mags_err['r3']**2+mags_err['r3']**2), np.sqrt(mags_err['mrk1018']**2+mags_err['r3']**2)), names = ('r3', 'mrk1018'))
    
    #---------------------------- sigma AGN only errors---------------------------
    
    negs = stats['lower_1sig'] < 0 #find obs where the mag lower bound is negative
    pos = stats['lower_1sig'] > 0 #mag lower bound is positive = fine
    
    negmags = stats[negs] #separate negative lower count errors
    posmags = stats[pos]
    
    lolims = np.ones(len(negmags), dtype=bool) #create array for neg lower bounds
    
    #deal with negative count data points
    
    upper_lims = negmags['upper_3sig'] #use upper 3 sigma limits for neg mag points
    ulims_corr = -2.5*np.log10(upper_lims) - mags['r3'][negs] #convert upper limit to mags
    
    #deal with positive count points
    
    centre_corr = -2.5*np.log10(posmags['centre']) - mags['r3'][pos] #convert centre for positive agn only counts points to mags
    
    lower_mag_corr = -2.5*np.log10(posmags['lower_1sig'])-mags['r3'][pos] #correct upper/lower bounds for positive count points with zeropoint mag
    upper_mag_corr = -2.5*np.log10(posmags['upper_1sig'])-mags['r3'][pos]
    
    lower_error = centre_corr-lower_mag_corr #find errors (diffr between central value and bounds)
    upper_error = upper_mag_corr-centre_corr
    
    asymmetric_error = [lower_error,upper_error] #create asymmetric error array for non-neg points
    
    #--------------------------correct for zeropoint mag----------------------
    
    #plt.gca().invert_yaxis()
    plt.ylim(23.5, 16)
    plt.ylabel('u\'-band magnitude')
    plt.xlabel('decimal year')
    #plot combined mags
    plt.errorbar(obs['dec_year'], mags_corr['mrk1018_corr']+zp, yerr = mags_corerr['mrk1018'], fmt='.', color='purple', label = 'Mrk 1018')
    #plot positive agn only points
    plt.errorbar(posmags['dec_yr'], centre_corr+zp, yerr=asymmetric_error, fmt='.b', label = 'AGN only')\
    #plot upper limits for negative count points
    plt.errorbar(negmags['dec_yr'], ulims_corr+zp, yerr=0.3, lolims = lolims, fmt='_b')
    plt.title('STELLA lightcurve')
    plt.ticklabel_format(useOffset=False)
    plt.legend(loc = 'lower right')
    plt.savefig('lightcurve_1sig')
    plt.close()
    #plt.show()
    
    #save data points
    p_dps = Table(data = (posmags['dec_yr'], centre_corr+zp, asymmetric_error[0], asymmetric_error[1]), names = ('dec_yr', 'agn_abs', 'lower_err', 'upper_err'))
    ascii.write(p_dps, 'data_points_pos.csv', format = 'csv', overwrite = True)
    
    n_dps = Table(data = (negmags['dec_yr'], ulims_corr+zp, np.zeros(len(ulims_corr)), np.zeros(len(ulims_corr))), names = ('dec_yr', 'agn_abs', 'lower_err', 'upper_err'))
    ascii.write(n_dps, 'data_points_neg.csv', format = 'csv', overwrite = True)
    
    dps = vstack([p_dps, n_dps])
    ascii.write(dps, 'data_points_all.csv', format = 'csv', overwrite = True)
