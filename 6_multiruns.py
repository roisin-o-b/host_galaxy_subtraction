import multiprocessing as mp
import glob
from astropy.io import ascii, fits
import numpy as np
import stella_pipeline_star

stella = ascii.read('/work1/roisin/work/mrk1018/optical/stella_pipeline/development+results/counts/counts_100621.csv')
vimos = fits.open('reproject/vimos_hostcut_reproj_20170730.fits')
bkgs = sorted(glob.glob('/work1/roisin/work/mrk1018/optical/stella_pipeline/development+results/local_bkgs/*.fits'))
ref3 = ascii.read('/work1/roisin/work/mrk1018/optical/stella_pipeline/development+results/host_subtraction/vimos_r3_counts.csv')
year = 2017.57939767283

i = tuple(np.arange(1000))

def multiple_runs(x):
    obs = stella_pipeline_star.host_sub(stella, bkgs, vimos, ref3, year)
    ascii.write(obs, 'runs/run'+str(x)+'.csv', format = 'csv', overwrite=True)

def pool_handler():
    p = mp.Pool()
    p.map(multiple_runs, i)
    p.close()
    p.join()
    
if __name__ == '__main__':
    
    pool_handler()
