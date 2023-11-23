from astropy.io import fits

vimos = fits.open('vimos_host_bkg_20170730.fits')[0]
vimos = vimos.data
data = fits.open('vimos_host_bkg_reproj_20170730.fits')[0]
data = data.data

from photutils.centroids import centroid_sources, centroid_com

i, j = centroid_sources(vimos, 201, 197, box_size=3, centroid_func=centroid_com)
x, y = centroid_sources(data, 1651, 1729, box_size=4, centroid_func=centroid_com)
data = data[round(y[0])-100:round(y[0])+100, round(x[0])-100:round(x[0])+100]
x, y = centroid_sources(data, 100, 100, box_size=4, centroid_func=centroid_com)

fits.writeto('vimos_hostcut_reproj_20170730.fits', data, overwrite = True)

#import sys
#sys.exit()

from photutils import CircularAperture, aperture_photometry

p1 = (i[0], j[0])
a1 = CircularAperture(p1, r = 10/0.205)

p2 = (x[0], y[0])
a2 = CircularAperture(p2, r = 10/0.322)

q = aperture_photometry(vimos, a1)
w = aperture_photometry(data, a2)

print(q, w)

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

norm = simple_norm(data, 'sqrt', percent=99.9)
plt.imshow(vimos, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
a1.plot(color='red', lw=1.5)
plt.xlim(0, vimos.shape[1]-1)
plt.ylim(0, vimos.shape[0]-1)
plt.show()

#data = data[round(y[0])-100:round(y[0])+100, round(x[0])-100:round(x[0])+100]
plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
a2.plot(color='red', lw=1.5)
plt.xlim(0, data.shape[1]-1)
plt.ylim(0, data.shape[0]-1)
plt.show()
