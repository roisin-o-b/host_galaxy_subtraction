import glob
import numpy as np
from astropy.io import ascii
import lc_fns as lc


files = sorted(glob.glob('runs/run*.csv')) #read in 1000 runs of stella pipeline

lc.group_dps(files) #call function to group agn_only count values by decimal year

######################################################################################

dpts = glob.glob('data_points/*.csv') #read in agn only counts, 1000 runs for each obs
dec_yrs = np.loadtxt('dec_yrs.csv') #load decimal year values

lc.agnonlystats(dpts, dec_yrs) #call function to retreive central and 1 sigma values from each obs

######################################################################################

obs = ascii.read('/home/roisin/work/mrk1018/optical_monitoring/stella_data/counts/obs_100621.csv') #read in STELLA data
stats = ascii.read('agnonlystats.csv', format='csv') #read in stat data from previous function

lc.plot_agnstats(obs, stats) #plot lightcurve and save values
