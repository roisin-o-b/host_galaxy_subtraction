import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, vstack
from astropy.convolution import convolve_fft
from astropy.visualization import simple_norm, SqrtStretch
from astropy.stats import sigma_clipped_stats
from photutils.centroids import centroid_sources, centroid_com
from photutils import detect_threshold, detect_sources, CircularAperture, aperture_photometry, Background2D

def host_sub(counts, loc_bkgs, host, host_star, host_date):
    
    host = host[0].data

    #create table for stella scaled host counts
    host_counts = Table(names = ('dec_yr', 'host_count'))
    
    for f, value in enumerate(counts):
        
        #extract stella seeing and convert to pixels
        fwhm_pix = counts['fwhm'][f]/0.322 #STELLA pixel scale
        
        #create gaussian based on stella seeing fwhm
        x, y = np.meshgrid(np.linspace(-30,30,60), np.linspace(-30,30,60))
        d = np.sqrt(x*x+y*y)
        sigma = fwhm_pix/(2*np.sqrt(2*np.log(2)))
        mu = 0.0
        g = np.exp( -( (d-mu)**2 / ( 2.0 * sigma**2) ) )
        
        # convolve with host image
        host_conv = convolve_fft(host, g, normalize_kernel = True)
        
        host_bkg = host_conv
        
        #read in local STELLA background
        loc_bkg = fits.open(loc_bkgs[f])
        loc_bkg = loc_bkg[0].data
        
        host_bkg = host_conv + loc_bkg
        
        #add random fluctuations to convolved host image (vimos + stella bkg)
        rows = 200
        columns = 200
        for j in range(columns):
            for i in range(rows):
                count = host_conv[i,j]
                del_c = np.sqrt(count)
                stbkg = loc_bkg[i,j]
                del_b = np.sqrt(stbkg)
                
                combined = np.sqrt(del_b**2 + del_c**2)
                
                pix_noise = np.random.default_rng().normal(0, combined, 1)
                host_bkg[i,j] = count + stbkg + pix_noise
        
        
        #subtract the background from the host
        host_map = Background2D(host_bkg, box_size=(40, 35))
        host_bkgsub = host_bkg - host_map.background
        
        #norm = simple_norm(host_bkgsub, percent=99)
        #plt.imshow(host_bkgsub, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
        #plt.xlim(0, host_bkgsub.shape[1]-1)
        #plt.ylim(0, host_bkgsub.shape[0]-1)
        #plt.savefig('host_bkgsub')
        #plt.close()
        
        #find host image centroid
        i, j = centroid_sources(host_bkgsub, 100, 100, box_size=5, centroid_func=centroid_com)
        
        #specify position and aperture
        position = (i[0], j[0])
        aperture = CircularAperture(position, r = 10/0.322)
        
        #do photometry
        host_tbl = aperture_photometry(host_bkgsub, aperture)
        host_counts.add_row([counts['dec_year'][f], host_tbl['aperture_sum']])
        
        #define vimos ref star 3 counts
        star_bool = host_star['dec_year'] == host_date
        star = host_star[star_bool]['ref3']

    ratio = counts['ref3'] / star
    ratio_er = ratio * np.sqrt( (counts['err3']/counts['ref3'])**2 + ((np.sqrt(star)/star)**2 ))

    host_pher = np.sqrt(host_counts['host_count'])

    counts['host'] = host_counts['host_count']*ratio #scale host according to ref star ratio

    counts['host_er'] = counts['host'] * np.sqrt( (host_pher/host_counts['host_count'])**2 + (ratio_er/ratio)**2 ) #propagate errors

    counts['agn_only'] = counts['mrk1018'] - counts['host']
    counts['err_agnonly'] = np.sqrt((counts['err_mrk1018'])**2 + counts['host_er']**2)
    
    #create table for averaged stella counts
    obs = Table(names = ('dec_yr', 'av_host', 'errhost', 'av_agnonly', 'err_agnonly'))

    #calculate error-weighted averages
    dates = np.unique(counts['dec_year'])

    for d in dates:
        bln = counts['dec_year'] == d #boolean for all exposures on the same day
        avH = np.average(counts[bln]['host'],       weights= 1/counts[bln]['host_er']**2)
        avAGN = np.average(counts[bln]['agn_only'], weights= 1/counts[bln]['err_agnonly']**2)
        
        errH = np.sqrt(1/np.sum(1/counts[bln]['host_er']**2))
        errA = np.sqrt(1/np.sum(1/counts[bln]['err_agnonly']**2))
    
        obs.add_row([d, avH, errH, avAGN, errA])
        
    return obs
