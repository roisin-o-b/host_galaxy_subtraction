from astropy.io import fits

stella = fits.open('new_science20210131A-0253EXP0004.fits')[0]
head = fits.open('VIMOS_Mrk1018_U_2017-07-30.fits')[0]
head = head.header

vimos = fits.open('vimos_host_bkg_20170730.fits')[0]
vimos.header = head

from astropy.wcs import WCS
import matplotlib.pyplot as plt

from reproject import reproject_interp
array, footprint = reproject_interp(vimos, stella.header)
array = array * (0.322**2 / 0.205**2)

fits.writeto('vimos_host_bkg_reproj_20170730.fits', array, stella.header, overwrite=True)
